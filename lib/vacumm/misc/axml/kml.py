# -*- coding: utf8 -*-
"""Interface to KML"""

#from base import *
#import os, tempfile
#
#__all__ = ['kml', 'kmz', 'Point', 'Placemark']
#
#class kml(Base):
#   """Top class for all kml documents. All elements must be child of a kml instance."""
#   pure_string = True # Strings are not protected (like <title>This is the title</title>, not <title>"This...)
#   def __init__(self, *args, **kwargs):
#       kwargs['xmlns'] = 'http://earth.google.com/kml/2.1'
#       base.__init__(self, *args, **kwargs)
#
#class kmz(kml):
#   self._files = []
#   def __init__(self, *args, **kwargs):
#       # kmz still has the kml tag name
#       kml.__init__(self, tag='kml', *args, **kwargs)
#   def add_files(self, *files):
#       """Add files to the zip"""
#       for file in files:
#           if os.path.exists(file) and file not in self._files:
#               self._files.append(file)
#   def write(self, file):
#       # Create a temporary directory
#       tmpdir = tempfile.mkdtemp(prefix='python-kmz-')
#
#
#class Folder(NameDesc, base):
#   pass
#
#class NameDesc:
#   exec def_setget('name', 'Name')
#   exec def_setget('description', 'Description', default='')
#
#
#class Point(Base):
#   def coordinates_formatter(self, value, reverse=False):
#       if not reverse:
#           return '%f, %f, %f'%value
#       return eval('(%s)'%value)
#   def set_coordinates(self, lon, lat, dep=0):
#       """Set coordinates"""
#       self._set_("coordinates", (lon, lat, dep))
#   def get_coordinates(self):
#       """Get coordinates
#       @return: (lon, lat, dep)
#       """
#       return self._get_("coordinates")
#
#
#class Placemark(NameDesc, Base):
#   exec def_setget('name', 'Name')
#   exec def_setget('description', 'Description', default='')
#   def set_point(self, lon, lat, dep=0):
#       for node in self:
#           if node.tag == 'Point':
#               del node
#       self.append(Point((lon, lat, dep)))
#
#class Icon(Base):
#   exec def_setget('href', 'Url of the image')
#
#class LatLonBox(Base):
#   exec def_setget('north', 'North bound')
#   exec def_setget('south', 'South bound')
#   exec def_setget('east', 'East bound')
#   exec def_setget('west', 'West bound')
#   exec def_setget('rotation', 'Rotation', default=0.)
#
#
