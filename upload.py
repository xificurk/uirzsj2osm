#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Import of data prepared from UIR-ZSJ.

"""

from __future__ import print_function

__author__ = "Petr Morávek (xificurk@gmail.com)"
__copyright__ = "Copyright (C) 2009-2011 Petr Morávek"
__license__ = "GPL 3.0"

__version__ = "1.1"

import logging
from optparse import IndentedHelpFormatter, OptionGroup, OptionParser
import os.path
import sys

sys.path.insert(0, os.path.join(sys.path[0], "libs"))

from colterm import *
import uirzsj2osm
import uir

init_translation(os.path.join(os.path.dirname(__file__), "locale"), languages=["cs"])
uir.init_translation(os.path.join(os.path.dirname(__file__), "locale"), languages=["cs"])


if __name__ == "__main__":
    # Setup console output logging
    colored_handler = ColoredStreamHandler(fmt="%(asctime)s %(levelname)-8s %(name)s >> %(message)s", datefmt="%Y-%m-%d %X")
    rootlog = logging.getLogger("")
    rootlog.addHandler(colored_handler)
    rootlog.setLevel(logging.WARN)

    # Parse command line arguements
    optp = OptionParser(formatter=IndentedHelpFormatter(max_help_position=20), conflict_handler="resolve", usage="pouziti: %prog [options] osm_id|ref [osm_id|ref [...]]", version="%prog "+__version__)
    optp.add_option("-n", "--no-color", help="vypnout pouzivani barevneho vystupu", dest="color", action="store_false", default=True)
    optp.add_option("-q", "--quiet", help="logging level ERROR", dest="loglevel", action="store_const", const=logging.ERROR, default=logging.WARN)
    optp.add_option("-v", "--verbose", help="logging level INFO", dest="loglevel", action="store_const", const=logging.INFO)
    optp.add_option("-d", "--debug", help="logging level DEBUG", dest="loglevel", action="store_const", const=logging.DEBUG)
    optp.add_option("-D", "--Debug", help="logging level ALL", dest="loglevel", action="store_const", const=0)

    opts,args = optp.parse_args()
    rootlog.setLevel(opts.loglevel)
    use_color(opts.color)

    if len(args) == 0:
        rootlog.critical("Je treba zadat alespon jedno OSM ID relace, nebo LAU2 kod obce.")
        raise SystemExit(1)
    filename = os.path.join(os.path.dirname(__file__), "credentials.txt")
    if not os.path.isfile(filename):
        rootlog.critical("Soubor s prihlasovacimi udaji pro OSM nenalezen. Vytvor soubor credentials.txt, na prvni radek zadej svoje uzivatelske jmeno a na druhy heslo.")
        raise SystemExit(1)
    with open(filename) as fp:
        uirzsj2osm.api.username = fp.readline().rstrip("\r\n")
        uirzsj2osm.api.password = fp.readline().rstrip("\r\n")
    for ref in args:
        imp = uirzsj2osm.Import(ref)
        imp.commit()
