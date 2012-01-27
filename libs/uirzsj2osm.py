# -*- coding: utf-8 -*-
"""
Scripts for import of UIR-ZSJ to OSM.

"""

from __future__ import print_function

__author__ = "Petr Morávek (xificurk@gmail.com)"
__copyright__ = "Copyright (C) 2011 Petr Morávek"
__license__ = "LGPL 3.0"

__version__ = "0.5.0"

from collections import defaultdict
import logging
import os.path
from pyproj import Proj, transform
import re
from time import sleep
import unicodedata

from colterm import *
import uir
import osmapis

colors["create"] = ANSI.color("G", False, "")
colors["modify"] = ANSI.color("RG", False, "")
colors["delete"] = ANSI.color("R", False, "")
colors["notice"] = ANSI.color("G", True, "")
logging.addLevelName(31, "NOTICE")

# Logging
log = logging.getLogger("uirzsj2osm")
log.notice = lambda msg: log.log(31, msg)


# Geometry
jtsk = Proj("+proj=krovak +ellps=bessel +towgs84=570.8,85.7,462.8,4.998,1.587,5.261,3.56")
wgs84 = Proj("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

def jtsk2wgs(x, y):
    return transform(jtsk, wgs84, x, y)

def wgs2jtsk(lon, lat):
    return transform(wgs84, jtsk, lon, lat)


# OSM
qhs = osmapis.QHS()
oapi = osmapis.OverpassAPI()
api = osmapis.API(auto_changeset={"enabled": True})

class PlaceNode(osmapis.Node):
    @classmethod
    def from_jtsk(cls, xy, tags):
        attribs = {}
        attribs["lon"], attribs["lat"] = jtsk2wgs(*xy)
        return cls(attribs, tags)

    @property
    def xy(self):
        return wgs2jtsk(self.lon, self.lat)

    @property
    def name(self):
        return self.tags.get("name:cs") or self.tags.get("name")


def to_ascii(msg):
    return unicodedata.normalize('NFKD', unicode(msg)).encode('ascii', 'ignore').lower()


class Reporter(object):
    place_msg = {"OBCE": u"obce", "COBE": u"části obce", "ZSJ": u"ZSJ"}

    def __init__(self, uir):
        self.uir = uir

    def set_obce(self, osm_id, ref):
        self.osm_id = osm_id
        self.ref = ref
        log.notice(u"Pripravuji data pro {} (relace {}).".format(self.uir.data["OBCE"][ref]["name"], osm_id))

    def admin_centre_name_mismatch(self, admin_centre):
        log.error(u"OSM data uvadi jako admin_centre uzel {} ({}), ale UIR-ZSJ data uvádí {}.".format(admin_centre.id, admin_centre.name, self.uir.data["OBCE"][self.ref]["name"]))

    def add_place(self, place, match, entity=None, attr={}):
        if place is None:
            return None
        msg = u"Uzel {} {} ({}): {{}}.".format(self.place_msg[entity], place.tags["name"], place.tags["place"])
        if match is None:
            log.warn(msg.format(u"nenalezen v OSM datech"))
            return None
        else:
            parts = [u"odpovídá uzlu {} ({}, vzdálenost {}m)".format(match.id, match.tags["place"], attr["dist"])]
            if len(match.refs["way"]) > 0:
                parts.append(u"uzel {} je součástí cest {}".format(match.id, ", ".join(str(ref.id) for ref in match.refs["way"])))
            if len(match.refs["relation"]) > 0:
                parts.append(u"uzel {} je součástí relací {}".format(match.id, ", ".join(str(ref.id) for ref in match.refs["relation"])))
            if len(parts) <= 1 and len(attr["tags_odbl"]) == 0:
                if match.odbl_problems is not None:
                    parts.append(u"uzel {} není kompatibilní s ODbL (ale neobsahuje důležité tagy)".format(match.id))
                log.info(msg.format(u", ".join(parts)))
                return None
            else:
                if match.odbl_problems is not None:
                    parts.append(u"uzel {} není kompatibilní s ODbL".format(match.id))
                log.error(msg.format(u", ".join(parts)))
                return u", ".join(parts)

    def add_unmatched(self, places):
        unmatched = defaultdict(list)
        for node in places:
            self.add_place(None, node)
            unmatched[node.tags["place"]].append(u"{} ({})".format(node.id, node.tags["name"]))
        for place, items in unmatched.items():
            log.error(u"Uzly s tagem place={}, které nebyly nalezeny v UIR-ZSJ datech: {}.".format(place, ", ".join(items)))


class Places(osmapis.OSM):
    tags_harmless = re.compile("^(created_by|name|place|population|fixme|note|ref|lau[12]|(source|is_in).*)$", re.I)
    tags_ignore = re.compile("^(created_by|name|place|population|source(:population)?|ref:(cobe|zsj)|lau[12])$", re.I)

    dist_shift = 30
    dist_obce = 3000
    dist_cobe = 2000
    dist_zsj = 1000

    def export(self, container, data):
        self.odbl_problems()
        self.build_refs(container)
        self.export_obce(container, data["OBCE"])
        self.export_cobe(container, data["COBE"])
        self.export_zsj(container, data["ZSJ"])
        self.reporter.add_unmatched(self.nodes.values())
        return container | self

    def odbl_problems(self, wait=1):
        log.debug("Analyzuji licenci pro uzly sidel.")
        try:
            qhs.problems(self)
        except osmapis.APIError:
            log.warn("Nepodarilo se overit licencni problemy uzlu... zkusim znovu za {}s.".format(wait))
            sleep(wait)
            self.odbl_problems(min(300, wait*2))

    def build_refs(self, container):
        for node in self.nodes.values():
            node.refs = {"admin": [], "relation": [], "way": []}
        for relation in container.relations.values():
            for member in relation.members:
                ref = member["ref"]
                if ref in self.nodes and member["type"] == "node":
                    if member["role"] == "admin_centre" and relation.tags.get("type") == "boundary" and relation.tags.get("boundary") == "administrative":
                        self.nodes[ref].refs["admin"].append(relation)
                    else:
                        self.nodes[ref].refs["relation"].append(relation)
        for way in container.ways.values():
            for ref in way.nds:
                if ref in self.nodes:
                    self.nodes[ref].refs["way"].append(way)

    def export_obce(self, container, place):
        match = self.find_obce(place)
        if place.pop("shift", False):
            place["xy"] = (place["xy"][0], place["xy"][1] + self.dist_shift)
        self.export_place(container, place, match, "OBCE")

    def export_cobe(self, container, places):
        for place in places:
            match = self.find_cobe(place)
            if place.pop("shift", False):
                place["xy"] = (place["xy"][0], place["xy"][1] - self.dist_shift)
            self.export_place(container, place, match, "COBE")

    def export_zsj(self, container, places):
        for place in places:
            match = self.find_zsj(place)
            self.export_place(container, place, match, "ZSJ")

    def find_obce(self, place):
        if to_ascii(place["name"]) == to_ascii(self.admin_centre.name):
            return self.admin_centre
        self.reporter.admin_centre_name_mismatch(self.admin_centre)
        return self.find_by_name(place, self.dist_obce) or self.find_by_ascii_name(place, self.dist_obce)

    def find_cobe(self, place):
        return self.find_by_ref(place, self.dist_cobe, "ref:cobe") or self.find_by_name(place, self.dist_cobe) or self.find_by_ascii_name(place, self.dist_cobe)

    def find_zsj(self, place):
        return self.find_by_ref(place, self.dist_cobe, "ref:zsj") or self.find_by_name(place, self.dist_cobe) or self.find_by_ascii_name(place, self.dist_cobe)

    def find_by_ref(self, place, score, ref_key):
        place_ref = place[ref_key]
        match = None
        for node in self.nodes.values():
            if node.tags.get(ref_key) != place_ref:
                continue
            dist = uir.distance(node.xy, place["xy"])
            if dist > score:
                continue
            score = dist
            match = node
        return match

    def find_by_name(self, place, score):
        place_name = place["name"].lower()
        match = None
        for node in self.nodes.values():
            if node.name.lower() != place_name:
                continue
            dist = uir.distance(node.xy, place["xy"])
            if dist > score:
                continue
            score = dist
            match = node
        return match

    def find_by_ascii_name(self, place, score):
        place_name = to_ascii(place["name"])
        match = None
        for node in self.nodes.values():
            if to_ascii(node.name) != place_name:
                continue
            dist = uir.distance(node.xy, place["xy"])
            if dist > score:
                continue
            score = dist
            match = node
        return match

    def export_place(self, export, place, match, entity):
        tags = dict(place)
        tags.pop("COBE", None)
        tags.pop("radius", None)
        tags.pop("population")
        tags["source"] = "csu:uir-zsj"
        place = PlaceNode.from_jtsk(tags.pop("xy"), tags)
        if match is None:
            self.reporter.add_place(place, match, entity)
            export.add(place)
        else:
            self.remove(match)
            attr = {"tags_odbl": set(), "dist": int(round(uir.distance(place.xy, match.xy)))}
            if match.odbl_problems is not None:
                attr["tags_odbl"] = set(filter(lambda k: self.tags_harmless.match(k) is None, match.tags.keys()))
                if len(match.refs["way"]) == 0 and len(match.refs["relation"]) == 0 and len(attr["tags_odbl"]) == 0:
                    self.reporter.add_place(place, match, entity, attr)
                    export.add(place)
            if place not in export:
                if len(match.refs["way"]) > 0 or len(match.refs["relation"]) > 0 or match.odbl_problems is not None:
                    place.tags["fixme"] = self.reporter.add_place(place, match, entity, attr)
                    if match.odbl_problems is None:
                        self.merge_tags(match, place)
                    export.add(match)
                    export.add(place)
                else:
                    self.reporter.add_place(place, match, entity, attr)
                    self.merge_tags(match, place)
                    for k, v in match.attribs.items():
                        if k not in ("lat", "lon"):
                            place.attribs[k] = v
                    export.add(place)

            for adm_rel in match.refs["admin"]:
                for member in adm_rel.members:
                    if member["role"] == "admin_centre" and member["ref"] == match.id:
                        member["ref"] = place.id

        if entity == "OBCE" and match != self.admin_centre:
            for member in export.relations[self.rel_id].members:
                if member["type"] == "node" and member["role"] == "admin_centre":
                    member["ref"] = place.id
            self.admin_centre.refs["admin"].remove(export.relations[self.rel_id])

    def merge_tags(self, match, place):
        if place.tags["place"] in ("suburb", "neighbourhood") and match.tags["place"] in ("village", "hamlet", "isolated_dwelling"):
            place.tags["place"] = match.tags["place"]
        if match.tags.get("name:cs", "").lower() == place.tags["name"].lower():
            place.tags["name"] = match.tags["name"]
        for k, v in match.tags.items():
            if self.tags_ignore.match(k) is None:
                if k == "note" and "note" in place.tags:
                    place.tags["note"] = u"{}, {}".format(place.tags["note"], v)
                else:
                    place.tags[k] = v


osmapis.wrappers["node"] = PlaceNode


# UIR-ZSJ import
class Import(object):
    query_places = '<query type="node" into="places"><area-query ref="{}"/><has-kv k="place"/><has-kv k="name"/></query>'
    query_rels = '<id-query type="relation" ref="{}"/><recurse type="node-relation" from="places"/>'
    query_print = '<print mode="meta" order="quadtile"/>'
    query_simple = '<union>' + query_places + query_rels + ' <recurse type="node-way" from="places"/><recurse type="way-node"/></union>' + query_print
    query_full = '<union>' + query_places + ' <union into="rels">' + query_rels + '</union><union><recurse type="node-way" from="places"/><recurse type="relation-way" from="rels"/></union><recurse type="way-node"/><recurse type="relation-node" from="rels"/></union>' + query_print

    def __init__(self, ref):
        if ref.isdigit():
            relation = self.download_by_id(int(ref))
        elif ref.startswith("CZ0"):
            relation = self.download_by_ref(ref)
        else:
            log.critical("{} neni ani OSM ID relace, ani hodnota ref tagu obce.")
            raise SystemExit(1)
        self.rel_id = relation.id
        self.rel_ref = relation.tags["ref"]
        self.uir = uir.Places(self.rel_ref)
        self.reporter = Reporter(self.uir)
        Places.reporter = self.reporter

    def download_by_ref(self, ref):
        query = '<query type="relation"><has-kv k="type" v="boundary"/><has-kv k="boundary" v="administrative"/><has-kv k="admin_level" v="8"/><has-kv k="ref" v="{}"/></query><print order="quadtile"/>'.format(ref)
        try:
            osm = oapi.interpreter(query)
        except osmapis.APIError as e:
            log.exception(e)
            log.critical("Nepodarilo se stahnout data relace s ref={}.".format(ref))
            raise SystemExit(1)
        if len(osm) == 1:
            return osm.pop()
        elif len(osm) == 0:
            log.critical("Relace s ref={} nenalezena.".format(ref))
        elif len(osm) > 1:
            log.critical("Urceni relace ref={} neni jednoznacne... vyber jedno z OSM id: {}".format(ref, ", ".join((str(osm_id) for osm_id in osm.relations))))
        raise SystemExit(1)

    def download_by_id(self, osm_id):
        query = '<id-query type="relation" ref="{}"/><print order="quadtile"/>'.format(osm_id)
        try:
            osm = oapi.interpreter(query)
        except osmapis.APIError as e:
            log.exception(e)
            log.critical("Nepodarilo se stahnout data relace {}.".format(osm_id))
            raise SystemExit(1)
        if osm_id not in osm.relations:
            log.critical("Relace {} nenalezena.".format(osm_id))
            raise SystemExit(1)
        rel = osm.relations[osm_id]
        if not (rel.tags.get("type") == "boundary" and rel.tags.get("boundary") == "administrative" and rel.tags.get("admin_level") == "8" and rel.tags.get("ref", "").startswith("CZ0")):
            log.critical("Relace {} neni relaci hranice obce.".format(osm_id))
            raise SystemExit(1)
        return rel

    def download_osm_data(self, simple=False):
        log.debug("Stahuji data pro relaci {}.".format(self.rel_id))
        if simple:
            query = self.query_simple.format(3600000000 + self.rel_id, self.rel_id)
        else:
            query = self.query_full.format(3600000000 + self.rel_id, self.rel_id)
        data = oapi.interpreter(query)
        places = Places(node for node in data.nodes.values() if "place" in node.tags and "name" in node.tags)
        places.admin_centre = self.find_admin_centre(data)
        log.info("V OSM datech nalezeno {} pojmenovanych uzlu s tagem place.".format(len(places)))
        return data, places

    def find_admin_centre(self, data):
        admin_centre = None
        for member in data.relations[self.rel_id].members:
            if member["type"] == "node" and member["role"] == "admin_centre":
                admin_centre = member["ref"]
                break
        return data.nodes.get(admin_centre)

    def prepare(self, simple=False):
        if simple:
            log.warn("POZOR: Nestahuji uzly a cesty odkazovane relacemi, pro editaci dat NENI MOZNE pouzit Merkaartor.")
        data, places = self.download_osm_data(simple)
        data.save(os.path.join(os.path.dirname(__file__), "..", "tmp", "{}.osm".format(self.rel_id)))
        data -= places
        places.rel_id = self.rel_id
        places.reporter.set_obce(self.rel_id, self.rel_ref)
        export = places.export(osmapis.OSM(data), self.uir.get_places(self.rel_ref))
        filename = os.path.join(os.path.dirname(__file__), "..", "import", "{}.osm".format(self.rel_id))
        export.save(filename)
        return os.path.normpath(filename)

    def commit(self):
        self.original = osmapis.OSM.load(os.path.join(os.path.dirname(__file__), "..", "tmp", "{}.osm".format(self.rel_id)))
        export = osmapis.OSM.load(os.path.join(os.path.dirname(__file__), "..", "import", "{}.osm".format(self.rel_id)))
        self.add_population(export)
        log.notice(u"Odesilam data pro {} (relace {}).".format(self.uir.data["OBCE"][self.rel_ref]["name"], self.rel_id))
        osc = osmapis.OSC.from_diff(self.original, export)
        self.upload(osc)

    def add_population(self, data):
        places = osmapis.OSM(node for node in data.nodes.values() if node.tags.get("place") in ("city", "town", "village", "hamlet") and "name" in node.tags)
        admin_centre = self.find_admin_centre(data)
        if admin_centre is None:
            log.critical("Relaci {} chybi admin_centre.".format(self.rel_id))
            raise SystemExit(1)
        pop_reduce = defaultdict(int)

        for node in places.nodes.values():
            if node == admin_centre:
                continue
            if len(set(node.tags.get("ref:cobe", "").split(";")) & set(self.uir.data["COBE"][self.rel_ref])) > 0:
                continue
            if len(set(node.tags.get("ref:zsj", "").split(";")) & set(self.uir.data["ZSJ"][self.rel_ref])) > 0:
                continue
            if ("ref:cobe" in node.tags or "ref:zsj" in node.tags) and node.id in self.original.nodes and "population" in self.original.nodes[node.id].tags:
                node.tags["population"] = self.original.nodes[node.id].tags["population"]
            else:
                log.warn(u"Nepodarilo se doplnit population tag pro uzel {} ({}) place={}.".format(node.id, node.tags["name"], node.tags["place"]))
            places.remove(node)

        for node in places.nodes.values():
            if node == admin_centre:
                continue
            if len(set(node.tags.get("ref:cobe", "").split(";")) & set(self.uir.data["COBE"][self.rel_ref])) > 0:
                continue
            pop = 0
            for ref in node.tags["ref:zsj"].split(";"):
                uir_node = self.uir.data["ZSJ"][self.rel_ref][ref]
                pop += uir_node["population"]
                pop_reduce[uir_node["COBE"]] += uir_node["population"]
            pop_reduce["OBCE"] += pop
            node.tags["population"] = str(pop)
            log.info(u"Doplnuji population={} pro uzel ZSJ {} ({}).".format(node.tags["population"], node.id, node.tags["name"]))
            places.remove(node)

        for node in places.nodes.values():
            if node == admin_centre:
                continue
            pop = 0
            for ref in node.tags["ref:cobe"].split(";"):
                uir_node = self.uir.data["COBE"][self.rel_ref][ref]
                pop += uir_node["population"] - pop_reduce.get(ref, 0)
            if pop < 0:
                log.error(u"Doplneni population pro uzel casti obce {} ({}) selhalo kvuli zaporne hodnote {}, kde udelali soudruzi z NDR chybu?".format(node.id, node.tags["name"], pop))
            else:
                node.tags["population"] = str(pop)
                pop_reduce["OBCE"] += pop
                log.info(u"Doplnuji population={} pro uzel casti obce {} ({}).".format(node.tags["population"], node.id, node.tags["name"]))
            places.remove(node)

        pop = self.uir.data["OBCE"][self.rel_ref]["population"] - pop_reduce["OBCE"]
        if pop < 0:
            log.error(u"Doplneni population pro uzel obce {} ({}) selhalo kvuli zaporne hodnote {}, kde udelali soudruzi z NDR chybu?".format(admin_centre.id, admin_centre.tags["name"], pop))
        elif pop > 0:
            admin_centre.tags["population"] = str(pop)
            log.info(u"Doplnuji population={} pro uzel obce {} ({}).".format(admin_centre.tags["population"], admin_centre.id, admin_centre.tags["name"]))

    def upload(self, data):
        print("Prehled importu:", color="header")
        for action_msg, action in (("VYTVORIT", "create"), ("UPRAVIT", "modify"), ("SMAZAT", "delete")):
            for elem_msg, elem in (("relace", "relations"), ("cesty", "ways"), ("uzly", "nodes")):
                container = getattr(getattr(data, action), elem)
                if len(container) > 0:
                    print("  {} {}: ".format(action_msg, elem_msg), end="", color=action)
                    parts = []
                    for element in container.values():
                        if "name" in element.tags:
                            parts.append(u"{} ({})".format(element.id, element.tags["name"]))
                        else:
                            parts.append(str(element.id))
                    print(", ".join(parts))
        if not prompt("Pokracovat s uploadem?", default="n", validate="BOOLEAN"):
            return
        try:
            log.notice("Uploaduji data...")
            changeset = api.auto_changeset["tags"]["comment"] = u"Import sídel UIR-ZSJ pro obci {}.".format(self.uir.data["OBCE"][self.rel_ref]["name"])
            api.upload_diff(data)
        except osmapis.APIError as e:
            log.exception(e)
            log.critical("Nepodarilo se nahrat data do OSM.")
            raise SystemExit(1)
