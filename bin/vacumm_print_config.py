#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Print your configuration of VACUMM"""
print 'This script is deprecated. Please use "vacumm_config.py print" instead'

# Backward compatibility
import sys
sys.argv = [sys.argv[0].replace('print_config', 'config'), 'print']+sys.argv[1:]
import vacumm_config
