# -*- coding: utf8 -*-
"""Module to archive files (codes, images, etc)"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#

import os,re,zipfile,shutil
from tempfile import mkdtemp
from vacumm.misc import FileTree
import shutil

__all__ = ['Archive', 'FigureArchive', 'CodeArchive']

class Archive(FileTree):
    """Basic archive class, derived from a file tree"""

    def __init__(self,input_dir,**kwargs):

        kwargs.setdefault('scan',False)
        FileTree.__init__(self,input_dir,relative=True,**kwargs)

    def __call__(self,*args,**kwargs):
        self.zipto(*args,**kwargs)

    def zip_to(self,zip_file=None,quiet=False,append=False):
        """Create a zip with all files

        @keyparam zip_file: Zip file name or zipfile.ZipFile instance
        @keyparam append: Append to archive instead of erasing
        @keyparam quiet: Say nothing
        """
        if zip_file is None: zip_file = input_dir+'.zip'
        if not quiet: print 'Zipping directory %s to %s...'% (self._input_dir,zip_file)
        if not isinstance(zip_file,zipfile.ZipFile):
            mode = 'w'
            if os.path.exists(zip_file):
                if append:
                    mode = 'a'
                else:
                    os.remove(zip_file)
            self._zip = zipfile.ZipFile(zip_file,mode)
        else:
            self._zip = zip_file

        # Generic archiving
        self._archive_to_('zip',quiet)
        if not quiet: print 'Created zip archive', zip_file

    def copy_to(self,dstdir,quiet=False,append=False):
        archdir = os.path.join(dstdir,os.path.basename(self._input_dir))
        if os.path.exists(archdir):
            if not append:
                shutil.rmtree(archdir)
        else:
            print 'creating',archdir
            os.makedirs(archdir)
        self._archdir = archdir
        if not quiet: print 'Copying directory %s to %s...'% (self._input_dir,dstdir)
        self._archive_to_('copy',quiet)
        if not quiet: print 'Created archive copy', archdir


    def _archive_to_(self,archive_method,quiet,dstdir=None,workdir=None):

        # Temporary directory for file transformations
        targetdir = os.path.basename(self._input_dir)
        if archive_method == 'zip':
            archive_cmd = self._zip.write
            archdir = targetdir
        elif archive_method == 'copy':
            archive_cmd = shutil.copy
            archdir = self._archdir
        if workdir is None:
            workdir = mkdtemp(prefix='tmp-python-archive-')
        else:
            if not os.path.exists(workdir): os.makedirs(workdir)
        tmpdir = os.path.join(workdir,targetdir)
        if not os.path.exists(tmpdir): os.makedirs(tmpdir)

        # Loop on file with code filtering
        for relfile in self.file_list():

            infile = os.path.join(self._input_dir,relfile)
            self._transfile = os.path.join(tmpdir,relfile)
            archfile = os.path.join(archdir,relfile)
            checkfiles = [self._transfile]
            if archive_method == 'copy':
                checkfiles.append(archfile)
            for checkfile in checkfiles:
                lastdir = os.path.dirname(checkfile)
                if not os.path.exists(lastdir): os.makedirs(lastdir)

            # File transformation (does nothing by default)
            self.file_transform(infile)
            if self._transfile is None:
                self._transfile = infile

            # Archive file
            archive_cmd(self._transfile,archfile)
            if not quiet: print '  added',archfile

        shutil.rmtree(workdir)


    def file_transform(self,infile):
        """Abstract class for file transformations"""
        self._transfile = None


class FigureArchive(Archive):
    default_patterns = FileTree.default_patterns.copy()
    suffixes = ['jpg','jpeg','gif','png','bmp','tiff']
    default_patterns['patterns'] = []
    for suf in suffixes:
        suf = '\.%s$'%suf
        default_patterns['patterns'].extend([suf,suf.upper()])

class CodeArchive(Archive):
    pass
