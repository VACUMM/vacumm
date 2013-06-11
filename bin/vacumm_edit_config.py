#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Edit your user configuration file of VACUMM"""

# Arguments
from optparse import OptionParser
import sys, os, shutil
parser = OptionParser(#usage="Usage: %prog [options]", 
    description="Edit your user configuration file of VACUMM")
parser.add_option('-e', '--editor', action='store', dest='editor',
    help='alternative editor (defaults to environment variable EDITOR or VISUAL)')    
    
# Parse
(options, args) = parser.parse_args()

# Import
try:
    from vacumm.config import edit_user_conf_file
except ImportError, e:
    parser.error("Error importing vacumm:\n"+e.message)

# Print
edit_user_conf_file(editor=options.editor)
