# -*- coding: utf-8 -*-
"""
Scripts for loading data from UIR-ZSJ.

"""

from __future__ import print_function

__author__ = "Petr Morávek (xificurk@gmail.com)"
__copyright__ = "Copyright (C) 2011 Petr Morávek"
__license__ = "LGPL 3.0"

__version__ = "0.5.0"

from collections import defaultdict
from gettext import translation
import logging
from math import sqrt
import os.path
import re

from colterm import *
from dbfpy.dbf import Dbf


log = logging.getLogger('uir')
log.addHandler(logging.NullHandler())


############################################################
### Gettext                                              ###
############################################################

def init_translation(localedir=None, languages=None, fallback=True):
    """
    Initialize gettext translation.

    Arguments:
        localedir   --- Directory with locales (see gettext.translation).
        languages   --- List of languages (see gettext.translation).
        fallback    --- Return a NullTranslations (see gettext.translation).

    """
    global _, gettext, ngettext
    trans = translation("uir", codeset="utf-8", localedir=localedir, languages=languages, fallback=fallback)
    _ = gettext = trans.gettext
    ngettext = trans.ngettext

init_translation(localedir=os.path.join(os.path.dirname(__file__), "locale"))


############################################################
### UIR-ZSJ                                              ###
############################################################

def distance(p1, p2):
    return sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

class Places(object):
    """
    Places container (OBCE, COBE, ZSJ)
    """
    directory = os.path.join(os.path.dirname(__file__), "..", "uir-zsj")
    dist_match = 2
    datasets = {"COBE": "KOD_CAST", "ZSJ": "KOD_ZSJ"}

    columns = {}
    columns["name"] = {"OBCE": "NAZOB", "COBE": "NAZCOBE", "ZSJ": "NAZZSJ"}
    columns["population"] = {"OBCE": "OB01", "COBE": "OB01", "ZSJ": "OBZSJ01"}
    columns["xy"] = {"OBCE": "S{}_G", "COBE": "S{}_G", "ZSJ": "S{}ZSJ_G"}

    pop_city = 90000
    pop_town = 3000
    pop_village = 50
    pop_neighbourhood = 1
    fix_name = re.compile(u"(.*) \(([^)]+)\)", re.I)

    def __init__(self, ref=""):
        self.data = {}

        if ref != "":
            log.debug(_("Loading UIR-ZSJ data for {}.").format(ref))
        else:
            log.debug(_("Loading all UIR-ZSJ data."))

        self.load_dataset("OBCE", ref)
        for name in self.datasets:
            self.load_dataset(name, ref)

    def load_dataset(self, name, ref):
        col_name = self.columns["name"][name]
        col_pop = self.columns["population"][name]
        col_x = self.columns["xy"][name].format("X")
        col_y = self.columns["xy"][name].format("Y")

        ref_key = self.datasets.get(name)
        if ref_key is None:
            dataset = {}
        else:
            dataset = defaultdict(dict)

        if name == "ZSJ":
            filename = os.path.join(self.directory, "ZSJD.DBF")
        else:
            filename = os.path.join(self.directory, "{}.DBF".format(name))

        for row in Dbf(filename, True):
            row = row.asDict()
            place_ref = row["LAU1"] + row["ICOB"]
            if not place_ref.startswith(ref) or (name == "ZSJ" and row["DD"] != row["DIL"]):
                continue
            place = {"name": row[col_name].decode("cp852"), "population": row[col_pop], "xy": (row[col_x], row[col_y])}
            if name == "ZSJ":
                place["COBE"] = row["KOD_CAST"]
            match = self.fix_name.search(place["name"])
            if match is not None:
                place["name"] = match.group(1)
                place["note"] = match.group(2)
            if ref_key is None:
                dataset[place_ref] = place
            else:
                dataset[place_ref][row[ref_key]] = place
        if ref_key is None:
            log.info(_("Loaded {} nodes from UIR-ZSJ dataset {}.").format(len(dataset), name))
        else:
            log.info(_("Loaded {} nodes from UIR-ZSJ dataset {}.").format(sum(len(row) for row in dataset.values()), name))
        self.data[name] = dataset

    def get_radius(self, place):
        return 39 * place["population"]**0.39 + 200

    def get_places(self, ref):
        if ref not in self.data["OBCE"]:
            self.log.critical(_("This container does not have data for {}.").format(ref))
            raise SystemExit(1)

        place = dict(self.data["OBCE"][ref])
        place["radius"] = self.get_radius(place)
        if place["population"] >= self.pop_city:
            place["place"] = "city"
        elif place["population"] >= self.pop_town:
            place["place"] = "town"
        else:
            place["place"] = "village"

        places = {"OBCE": place}

        if "COBE" in self.data:
            places["COBE"] = []
            for place_ref, place in self.data["COBE"][ref].items():
                append = True
                is_separate = False
                if distance(place["xy"], places["OBCE"]["xy"]) < self.dist_match:
                    if place["name"].lower() == places["OBCE"]["name"].lower():
                        places["OBCE"]["ref:cobe"] = place_ref
                        append = False
                    else:
                        places["OBCE"]["shift"] = True
                else:
                    is_separate = distance(place["xy"], places["OBCE"]["xy"]) > places["OBCE"]["radius"]
                if append:
                    place = dict(place)
                    place["ref:cobe"] = place_ref
                    place["radius"] = self.get_radius(place)
                    if place["population"] == 0:
                        place["place"] = "locality"
                    elif is_separate:
                        if place["population"] > self.pop_village:
                            place["place"] = "village"
                        else:
                            place["place"] = "hamlet"
                    else:
                        if places["OBCE"]["place"] in ("city", "town"):
                            place["place"] = "suburb"
                        else:
                            place["place"] = "neighbourhood"
                    places["COBE"].append(place)

            if "ZSJ" in self.data:
                places["ZSJ"] = []
                for place_ref, place in self.data["ZSJ"][ref].items():
                    append = True
                    is_separate = False
                    if distance(place["xy"], places["OBCE"]["xy"]) < self.dist_match:
                        if place["name"].lower() == places["OBCE"]["name"].lower():
                            places["OBCE"]["ref:zsj"] = place_ref
                            append = False
                            places["OBCE"]["shift"] = False
                        else:
                            places["OBCE"]["shift"] = True
                    else:
                        is_separate = distance(place["xy"], places["OBCE"]["xy"]) > places["OBCE"]["radius"]
                    for node in places["COBE"]:
                        if distance(place["xy"], node["xy"]) < self.dist_match:
                            if place["name"].lower() == node["name"].lower():
                                node["ref:zsj"] = place_ref
                                append = False
                            else:
                                node["shift"] = True
                        if node["ref:cobe"] == place["COBE"] and is_separate:
                            is_separate = distance(place["xy"], node["xy"]) > node["radius"]
                    if append:
                        place = dict(place)
                        place["ref:zsj"] = place_ref
                        if place["population"] < self.pop_neighbourhood:
                            place["place"] = "locality"
                        elif is_separate:
                            if place["population"] > self.pop_village:
                                place["place"] = "village"
                            else:
                                place["place"] = "hamlet"
                        else:
                            place["place"] = "neighbourhood"
                        places["ZSJ"].append(place)

                log.info(_("Merging of nodes complete: UIR-ZSJ data contains 1 OBCE, {} COBE, {} ZSJ.").format(len(places["COBE"]), len(places["ZSJ"])))
            else:
                log.info(_("Merging of nodes complete: UIR-ZSJ data contains 1 OBCE, {} COBE.").format(len(places["COBE"])))

        return places
