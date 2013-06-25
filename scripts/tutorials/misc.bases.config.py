#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

from vacumm.misc.bases import Object

class MyObject(Object):
    # These are the defaults but in case you want to change these behavior...
    @classmethod
    def get_config_spec_file(cls):
        return os.path.splitext(__file__)[0] + '.ini'
        # => misc.bases.config.ini
    @classmethod
    def get_config_section_name(cls):
        # => MyObject
        return cls.__name__

def main():
    
    obj = MyObject()
    cfg = obj.get_default_config()
    obj.info('Default configuration:\n%s', obj.pformat(cfg.dict()))
    cfgfile = os.path.splitext(__file__)[0] + '.cfg'
    cfg = obj.load_config(cfgfile)
    obj.info('Loaded configuration (%s):\n%s', cfgfile, obj.pformat(cfg.dict()))
    obj.debug('Debug message after configuration has been loaded')

if __name__ == '__main__':
    main()

