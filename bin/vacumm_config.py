#!/usr/bin/env python
# -*- coding: utf8 -*-
"""VACUMM configuration"""

# Arguments
from argparse import ArgumentParser
import sys, os, shutil
try:
    from vacumm.config import handle_print_args, handle_edit_args
    from vacumm.__init__ import __version__
except ImportError, e:
    sys.exit("Error importing vacumm:\n"+e.message)

# Main parser
parser = ArgumentParser(#usage="Usage: %prog [options]",
    description="Print your configuration of VACUMM (%s)"%__version__)

# Sub parsers
subparsers = parser.add_subparsers()

# - print
parser_print = subparsers.add_parser('print', help='print configuration')
parser_print.add_argument('--system', '-s', action='store_true',
    help="Print only system information")
parser_print.add_argument('--direc', '-d', action='store_true',
    help="Print only important VACUMM directories")
parser_print.add_argument('--config', '-c', action='store_true',
    help="Print only VACUMM config itself")
parser_print.add_argument('--packages', '-p', action='store_true',
    help="Print only version of important packages")
parser_print.add_argument('--no-header', action='store_true',
    help="Do not add headers")
parser_print.set_defaults(func=handle_print_args)

# - edit
parser_edit = subparsers.add_parser('edit', help='edit configuration')
parser_edit.add_argument('-e', '--editor', action='store', dest='editor',
    help='alternative editor (defaults to environment variable EDITOR or VISUAL)')
parser_edit.set_defaults(func=handle_edit_args)


# Parse
args = parser.parse_args()
args.func(args)
