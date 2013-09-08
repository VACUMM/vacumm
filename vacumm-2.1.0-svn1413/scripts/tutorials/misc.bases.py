#!/usr/bin/env python
# -*- coding: utf-8 -*-

from vacumm.misc.bases import Object

class MyObject(Object):
    
    @classmethod
    def init_class(cls, name, bases, dct):
        # do your class related stuff here
        print 'stack trace:'
        print cls.stack_trace()
        print
        print 'default config:'
        print cls.get_default_config()
    
    def __init__(self):
        Object.__init__(self)
        # ...

def main():
    o = MyObject()

if __name__ == '__main__':
    main()

