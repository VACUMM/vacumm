#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from vcmq import Object


class MyObject(Object):

    @classmethod
    def init_class(cls, name, bases, dct):
        # do your class related stuff here
        print('stack trace:')
        print(cls.stack_trace())
        print()
        print('default config:')
        print(cls.get_default_config())

    def __init__(self):
        Object.__init__(self)
        # ...


MyObject()
