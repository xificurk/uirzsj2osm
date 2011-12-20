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


class Places(osmapis.OSM):
    tags_harmless = re.compile("^(created_by|name|place|population|lau[12]|(source|is_in).*)$", re.I)
    tags_ignore = re.compile("^(created_by|name|place|population|source(:population)?|ref:(cobe|zsj)|lau[12])$", re.I)

    dist_shift = 30
    dist_obce = 3000
    dist_cobe = 2000
    dist_zsj = 1000

    def export(self, container, data):
        self.stats = {"match": {}, "distance": defaultdict(list), "tags_odbl": defaultdict(int)}
        self.odbl_problems()
        self.build_refs(container)
        self.export_obce(container, data["OBCE"])
        self.export_cobe(container, data["COBE"])
        self.export_zsj(container, data["ZSJ"])
        unmatched = defaultdict(list)
        for node in self.nodes.values():
            self.stats["match"].setdefault(None, defaultdict(list))[node.tags["place"]].append(node.tags["name"])
            unmatched[node.tags["place"]].append(u"{} ({})".format(node.id, node.tags["name"]))
        for place, items in unmatched.items():
            log.error(u"Uzly s tagem place={}, ktere nebyly nalezeny v UIR-ZSJ datech: {}.".format(place, ", ".join(items)))
        return container | self

    def odbl_problems(self):
        log.debug("Analyzuji licenci pro uzly sidel.")
        qhs.problems(self)

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
        msg = u"Uzel obce {} ({}): {{}}.".format(place["name"], place["place"])
        self.export_place(container, place, match, msg)

    def export_cobe(self, container, places):
        for place in places:
            match = self.find_cobe(place)
            if place.pop("shift", False):
                place["xy"] = (place["xy"][0], place["xy"][1] - self.dist_shift)
            msg = u"Uzel casti obce {} ({}): {{}}.".format(place["name"], place["place"])
            self.export_place(container, place, match, msg)

    def export_zsj(self, container, places):
        for place in places:
            match = self.find_zsj(place)
            msg = u"Uzel ZSJ {} ({}): {{}}.".format(place["name"], place["place"])
            self.export_place(container, place, match, msg)

    def find_obce(self, place):
        if place["name"].lower() == self.admin_centre.tags["name"].lower() and uir.distance(place["xy"], self.admin_centre.xy) < self.dist_obce:
            return self.admin_centre
        log.error(u"OSM data uvadi jako admin_centre uzel {} ({}), ale UIR-ZSJ data rikaji {}.".format(self.admin_centre.id, self.admin_centre.tags["name"], place["name"]))
        score = self.dist_obce
        match = None
        for node in self.nodes.values():
            if place["name"].lower() != node.tags["name"].lower():
                continue
            dist = uir.distance(place["xy"], node.xy)
            if dist > score:
                continue
            score = dist
            match = node
        return match

    def find_cobe(self, place):
        return self.find_cobe_zsj(place, self.dist_cobe, "ref:cobe")

    def find_zsj(self, place):
        return self.find_cobe_zsj(place, self.dist_zsj, "ref:zsj")

    def find_cobe_zsj(self, place, score, ref_key):
        match = None
        for node in self.nodes.values():
            if node.tags.get(ref_key) != place[ref_key]:
                continue
            dist = uir.distance(node.xy, place["xy"])
            if dist > score:
                continue
            score = dist
            match = node
        if match is not None:
            return match
        for node in self.nodes.values():
            if node.tags.get("name").lower() != place["name"].lower():
                continue
            dist = uir.distance(node.xy, place["xy"])
            if dist > score:
                continue
            score = dist
            match = node
        return match

    def export_place(self, export, place, match, msg):
        tags = dict(place)
        tags.pop("COBE", None)
        tags.pop("radius", None)
        tags.pop("population")
        tags["source"] = "csu:uir-zsj"
        place = PlaceNode.from_jtsk(tags.pop("xy"), tags)
        if match is None:
            self.stats["match"].setdefault(place.tags["place"], defaultdict(list))[None].append(place.tags["name"])
            log.warn(msg.format("nenalezen v OSM datech"))
            export.add(place)
        else:
            self.remove(match)
            d = int(round(uir.distance(place.xy, match.xy)))
            self.stats["match"].setdefault(place.tags["place"], defaultdict(list))[match.tags["place"]].append(place.tags["name"])
            self.stats["distance"][place.tags["place"]].append(d)
            parts = ["odpovida uzlu {} ({}, vzdalenost {}m)".format(match.id, match.tags["place"], d)]
            if len(match.refs["way"]) > 0:
                parts.append("uzel {} je soucasti cest {}".format(match.id, ", ".join(str(ref) for ref in match.refs["way"])))
            if len(match.refs["relation"]) > 0:
                parts.append("uzel {} je soucasti relaci {}".format(match.id, ", ".join(str(ref) for ref in match.refs["relation"])))
            if match.odbl_problems is not None:
                relevant_tags = set(filter(lambda k: self.tags_harmless.match(k) is None, match.tags.keys()))
                if len(parts) <= 1 and len(relevant_tags) == 0:
                    export.add(place)
                    parts.append("uzel {} neni kompatibilni s ODbL (ale neobsahuje zadne dulezite tagy)".format(match.id))
                    log.info(msg.format(", ".join(parts)))
                else:
                    parts.append("uzel {} neni kompatibilni s ODbL".format(match.id))
                    if len(relevant_tags) > 0:
                        self.stats["tags_odbl"][None] += 1
                        for tag in relevant_tags:
                            self.stats["tags_odbl"][tag] += 1
            if place not in export:
                if len(parts) > 1:
                    if match.odbl_problems is None:
                        self.merge_tags(match, place)
                    place.tags["fixme"] = ", ".join(parts)
                    export.add(match)
                    export.add(place)
                    log.error(msg.format(", ".join(parts)))
                else:
                    self.merge_tags(match, place)
                    for k, v in match.attribs.items():
                        if k not in ("lat", "lon"):
                            place.attribs[k] = v
                    export.add(place)
                    log.info(msg.format(", ".join(parts)))

            for adm_rel in match.refs["admin"]:
                for member in adm_rel.members:
                    if member["role"] == "admin_centre" and member["ref"] == match.id:
                        member["ref"] = place.id

    def merge_tags(self, match, place):
        if place.tags["place"] in ("suburb", "neighbourhood") and match.tags["place"] in ("village", "hamlet"):
            place.tags["place"] = match.tags["place"]
        for k, v in match.tags.items():
            if self.tags_ignore.match(k) is None:
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
        log.notice(u"Pripravuji data pro {} (relace {}).".format(self.uir.data["OBCE"][self.rel_ref]["name"], self.rel_id))
        export = places.export(osmapis.OSM(data), self.uir.get_places(self.rel_ref))
        filename = os.path.join(os.path.dirname(__file__), "..", "import", "{}.osm".format(self.rel_id))
        export.save(filename)
        return os.path.normpath(filename)

    def commit(self):
        original = osmapis.OSM.load(os.path.join(os.path.dirname(__file__), "..", "tmp", "{}.osm".format(self.rel_id)))
        export = osmapis.OSM.load(os.path.join(os.path.dirname(__file__), "..", "import", "{}.osm".format(self.rel_id)))
        self.add_population(export)
        log.notice(u"Odesilam data pro {} (relace {}).".format(self.uir.data["OBCE"][self.rel_ref]["name"], self.rel_id))
        osc = osmapis.OSC.from_diff(original, export)
        self.upload(osc)

    def add_population(self, data):
        places = osmapis.OSM(node for node in data.nodes.values() if node.tags.get("place") in ("city", "village", "hamlet") and "name" in node.tags)
        admin_centre = self.find_admin_centre(data)
        if admin_centre is None:
            log.critical("Relaci {} chybi admin_centre.".format(self.rel_id))
            raise SystemExit(1)

        for node in places.nodes.values():
            if node == admin_centre or "ref:cobe" in node.tags or "ref:zsj" in node.tags:
                continue
            log.warn(u"Nepodarilo se doplnit population tag pro uzel {} ({}) place={}.".format(node.id, node.tags["name"], node.tags["place"]))
            places.remove(node)
        reduce_cobe = defaultdict(int)
        for node in places.nodes.values():
            if node == admin_centre or "ref:cobe" in node.tags:
                continue
            uir_node = self.uir.data["ZSJ"][self.rel_ref][node.tags["ref:zsj"]]
            node.tags["population"] = str(uir_node["population"])
            reduce_cobe[uir_node["COBE"]] += uir_node["population"]
            log.info(u"Doplnuji population={} pro uzel ZSJ {} ({}).".format(node.tags["population"], node.id, node.tags["name"]))
            places.remove(node)
        reduce_obce = sum(val for val in reduce_cobe.values())
        for node in places.nodes.values():
            if node == admin_centre:
                continue
            pop = self.uir.data["COBE"][self.rel_ref][node.tags["ref:cobe"]]["population"] - reduce_cobe.get(node.tags["ref:cobe"], 0)
            if pop < 0:
                log.error(u"Doplneni population pro uzel casti obce {} ({}) selhalo kvuli zaporne hodnote {}, kde udelali soudruzi z NDR chybu?".format(node.id, node.tags["name"], pop))
            else:
                node.tags["population"] = str(pop)
                reduce_obce += pop
                log.info(u"Doplnuji population={} pro uzel casti obce {} ({}).".format(node.tags["population"], node.id, node.tags["name"]))
            places.remove(node)
        pop = self.uir.data["OBCE"][self.rel_ref]["population"] - reduce_obce
        if pop < 0:
            log.error(u"Doplneni population pro uzel obce {} ({}) selhalo kvuli zaporne hodnote {}, kde udelali soudruzi z NDR chybu?".format(admin_centre.id, admin_centre.tags["name"], pop))
        else:
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
