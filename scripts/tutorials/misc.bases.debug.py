#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cdms2, numpy

from vacumm.misc.axes import create_time, create_lat, create_lon
from vacumm.misc.bases import Object

class MyObject(Object):
    def do_something(self):
        num, den = 0, 0
        try: num / den
        except: self.exception('Division by 0 failed !\n%s', self.exception_trace())
        self.info('The function emitting this message is: %s', self.func_name())
        self.info('The stack trace of the function emitting this message is:\n%s', self.stack_trace())
        t = create_time(['1900-01-01 00:00:00', '9999-12-31 23:59:59'])
        y = create_lat(range(-90, 90))
        x = create_lon(range(-180, 180))
        v = cdms2.createVariable(numpy.random.ranf((len(t), len(y), len(x))), axes=[t, y, x], id='data', long_name='random data')
        self.info('Time:\n%s', self.describe(t))
        self.info('Latitude:\n%s', self.describe(y))
        self.info('Longitude:\n%s', self.describe(x, stats=True))
        self.info('Variable:\n%s', self.describe(v, stats=True))

def main():
    obj = MyObject()
    obj.do_something()

if __name__ == '__main__':
    main()

