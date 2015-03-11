#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Edit your user configuration file of VACUMM"""
print 'This script is deprecated. Please use "vacumm_config.py edit" instead'

# Backward compatibility
import sys
sys.argv = [sys.argv[0].replace('edit_config', 'config'), 'edit']+sys.argv[1:]
import vacumm_config
