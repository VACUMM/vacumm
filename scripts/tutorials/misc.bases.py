#!/usr/bin/env python
# -*- coding: utf-8 -*-


from vacumm.misc.bases import Object


class MyObject(Object):
    
    @classmethod
    def init_class(cls, name, bases, dct):
        super(Object, cls).init_class(name, bases, dct)
        # do your class related stuff here


def main():
    pass


if __name__ == '__main__':
    main()


