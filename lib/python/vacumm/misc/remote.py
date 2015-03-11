# -*- coding: utf8 -*-
"""Working with remote files

.. warning:: Yoo need `paramiko <http://www.lag.net/paramiko/>`_ to use this module.
"""
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

import os, shutil, re
from stat import *
import glob
from urlparse import urlparse
from warnings import filterwarnings
import subprocess
filterwarnings('ignore', message='.*Not using mpz_powm_sec.*')

__all__ = ['SSHBank', 'sshbank', 'WorkFile', 'InputWorkFiles', 'OutputWorkFile']

class SSHBank(object):
    """Interface to handle a bank of SSH/SFTP connections

    :Example:

        >>> sshbank = SSHBank()
        >>> client = sshbank('camarpor')
        >>> ssh = client('ssh')
        >>> res = ssh.exec_command('ls')
        >>> ftp = client('ftp')
        >>> ftp.get('remote_file', 'local_file')
        >>> ssh2 = sshbank.ssh('myhost:4022')
        >>> ftp2 = sshbank.ftp('user@myhost')
    """
    def __init__(self):
        try:
            import paramiko
        except ImportError:
            raise('You need paramiko to use SSHBank')
        self.paramiko = paramiko
        self._bank = {}
    def __call__(self, host):

        # Check the bank
        if self._bank.has_key(host):
            return self._bank[host]

        # Get host, etc
        hostid = host
        #print 'hostid',hostid
        ssh = self.paramiko.SSHClient()
        ssh.load_system_host_keys()
        pp = urlparse('sftp://'+host)
        username = pp.username
        port = pp.port if pp.port is not None else 22
        host = pp.hostname

        # Connect
        ssh.set_missing_host_key_policy(self.paramiko.AutoAddPolicy())
        ssh.load_system_host_keys()
        if port is not None: port = int(port)
        ssh.connect(host, port=port, username=username)
        ftp = ssh.open_sftp()
        self._bank[hostid] = dict(ssh=ssh, ftp=ftp)
        return self._bank[hostid]

    def ssh(self, host):
        """Get the SSH agent"""
        return self(host)['ssh']

    def ftp(self, host):
        """Get the SFTP agent"""
        return self(host)['ftp']

try:
    #: Ready to use SSH bank
    sshbank = SSHBank()
except:
    sshbank = None


class WorkFileException(Exception):
    pass

class WorkFile(object):
    """Base class for :class:`InputWorkFiles` and class:`OutputWorkFile`

    :Params:

        - **logger**:  A :class:`~vacumm.misc.io.Logger` (or subclass) instance
        - **ssh**: A ssh connexion (for instance created using :class`SSHBank`
          connected to an host).
        - **umask**, optional: Argument to :func:`os.umask`
        - **dmode**, optional: Directory unix mode (see :func:`os.chmod`)
        - **fmode**, optional: File unix mode (see :func:`os.chmod`)
    """
    def __init__(self, logger, ssh, umask=0, dmode=S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH,
        fmode=S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH, raise_error=True):
        self._raise_error = raise_error
        self.logger  = logger
        if ssh is not False:
            self.ssh = ssh['ssh']
            self.ftp = ssh['ftp']
        else:
            self.ssh = self.ftp = None
        if umask is not None: os.umask(umask)
        self.fmode = fmode
        self.dmode = dmode

    def error(self, msg, warn=None):
        """Send an ERROR message and exit"""
        if warn is None: warn = not self._raise_error
        if self.logger:
            if warn:
                self.logger.warning(msg)
            else:
                self.logger.error(msg)
                raise WorkFileException(msg)
        else:
            if not warn:
                raise WorkFileException(msg)
            else:
                print msg
    def debug(self, msg):
        """Send a DEBUG message"""
        if self.logger:
            self.logger.debug(msg)
        else:
            print msg
    def warning(self, msg):
        """Send a WARNING message"""
        if self.logger:
            self.logger.warning(msg)
        else:
            print msg
    def info(self, msg):
        """Send an INFO message"""
        if self.logger:
            self.logger.info(msg)
        else:
            print msg

    def remote_exec(self, cmd):
        """Execute a remote command and return result as list of lines

        :Example:

            >>> files = workfile.remote_exec('ls')

        """
        if self.ssh is None:
            self.warning('No remote command in local mode')
        else:
            res = self.ssh.exec_command(cmd)
            return [line[:-1] for line in res[1].readlines()]

    def _host2ssh_(self, sshbank):
        if  self.host:
            if sshbank is None:
                from vacumm.misc.remote import sshbank
                #sshbank = SSHBank()
                if sshbank is None:
                    raise WorkFileException("Can't import sshbank because paramiko is not found")
            return sshbank(self.host)
        return False

    @classmethod
    def expand_path(cls, path, subst=None):
        """Expand '~', environment and other variables in a path

        :Params:

            - **path**: Path to be expanded
            - **subst**, optional: Dictionnary of variables that may be used
              for expansion
        """
        exps = globals().copy()
        exps.update(os.environ)
        for vn in "user", "home":
            exps.setdefault(vn, os.environ[vn.upper()])
        if subst is not None:
            exps.update(subst)
        return os.path.expanduser(path%exps)


    @classmethod
    def parse_path(cls, path, rootdir=None, subst=None, mode=None):
        """Get the (inputdir, outputdir, path) from path (and rootdir)

        :Params:

            - **path**: Simple or complex path with remote and local part.
              The generic form is :
              ``"[<prefix>][(<remote_dir>><local_dir>)[<pattern>]"``
              where ``<pattern>`` is global expression pattern (see :func:`glob.glob`).
              It can also take the form of a tuple.
              Remote paths follow specifications from RFC1738:
              http://tools.ietf.org/html/rfc1738.html.

              .. warning::

                  Be careful when using ``<prefix>`` with a remote path.

            - **rootdir**, optional: Optional prefix to prepend to path.
            - **subst**, optional: Dictionary of variables used for string isubstitutions on path.
              All environnement variables (taken from :attr:`os.environ`) are also
              substituted.

        :Example:

            >>> obj.parse_path('(sftp://user@host.fr:1022/my/path>/local/path)/to/data.nc', mode='get')
            >>> obj.parse_path(('sftp://user@host.fr:1022/my/path','/local/path','/to/data.nc'), mode='get')
            >>> obj.obj.parse_path('/local/path/data.nc')
            >>> obj.obj.parse_path('/home10(/user1>/user2)/data.nc')

        """
        # Root prefix
        if isinstance(path, tuple):
            if len(path)==0:
                path = ''
            elif len(path)==1:
                path = path[0]
            elif len(path)==2:
                path = '(%s,%s)'%path
            elif len(path)==3:
                path = '(%s>%s)%s'%path
            else:
                path = '%s(%s>%s)%s'%path[:4]
        else:
            path = path.strip()
        if rootdir is not None:
            path = os.path.join(rootdir.strip(), path)

        # Check '>' path separators
        if path.find('>') != -1 and path.rfind('>') != path.find('>'):
            raise WorkFileException('There is more than one ">" in the specified formatted path: '+path)

        # Substitutions
        path = cls.expand_path(path, subst=subst)

        # Split parts
        parsed = re_workdir_parse(path)
        if not parsed:
            return None,None,path.strip()
        spath = list(parsed[0])
        indir = os.path.join(spath[0], spath[1])
        outdir = os.path.join(spath[0], spath[2])
        #if '>' in indir:
            #raise WorkFileException('Input path must not contain a ">": %s'%indir)
        #if '>' in outdir:
            #raise WorkFileException('Output path must not contain a ">": %s'%outdir)

        # Check input and output dirs
        if mode is None:
            if hasattr(cls, 'get'): mode = 'get'
            elif hasattr(cls, 'put'): mode = 'put'
        ppi = urlparse(indir)
        ppo = urlparse(outdir)
        if (mode == 'get' and ppo.hostname is not None) or \
           (mode == 'put' and ppi.hostname is not None):
            indir,outdir = outdir,indir
        pattern = spath[3].strip(os.path.sep+' ')
        return indir, outdir, pattern
re_workdir_parse = re.compile('^(.*)\(([^>]+)>([^>]+)\)(.*)$').findall

class InputWorkFiles(WorkFile):
    """A class to deal with input remote files

    :Params:

        - *path*: path in the form ``"<prefix>(<remote_dir>><local_dir>)<pattern>"``.
          (see :func:`parse_path`)

    :Params:

        - *logger*: a :class:`~vacumm.misc.io.Logger` (or subclass) instance
        - *sshbank*: a :class`SSHBank` instance
        - *transfer*: automatically start the transfer after initialization

    :Example:

    >>> wfile = InputWorkFiles('(caparmor-sftp:/home125>/home200/caparmor)toto*/data/file*.nc')
    >>> wfile.get() # update
    >>> wfile = InputWorkFiles('/home510/toto/toto.nc') # does nothing !
    >>> wfile = InputWorkFiles(('/home15>/home12)toto/toto.nc') # local copy only
    >>> print wfile.local_files()
    >>> print wfile.remote_files()

    """
    def __init__(self, path, rootdir=None, logger=None, sshbank=None, transfer=False, subst=None, check=2, sort=None, filter=None, raise_error=True, **kwargs):
        # Parse path
        self._raise_error = raise_error
        indir,outdir,path = self.parse_path(path, rootdir=rootdir, subst=subst)
        #print 's3',indir,outdir,path
        self.sort = sort
        self.filter = filter
        if indir is None: # pure local with nothing to copy
            self.host = ''
            self.local_pattern = self.remote_pattern = path
            self.cpmode = 0
            self.rem2loc = self.loc2rem = None
            WorkFile.__init__(self, logger, False, raise_error=raise_error, **kwargs)
            return
        self.remote_dir, self.local_dir, self.pattern = indir,outdir,path
        self.local_pattern = os.path.join(self.local_dir, self.pattern)

        # Extract host and port
        up = urlparse(self.remote_dir)
        if up.hostname is None:
            self.host = ''
            self.cpmode = 1
        else:
            self.host = up.netloc
            self.cpmode = 2
        self.remote_dir = up.path
        self.remote_pattern = os.path.join(self.remote_dir, self.pattern)
        #self.rem2loc = {}
        #self.loc2rem = {}
        self._local_files = None

        # Get SSH agent from bank
        ssh = self._host2ssh_(sshbank)

        # Final init
        WorkFile.__init__(self, logger, ssh, raise_error=raise_error, **kwargs)

        # Already start the transfer?
        self.check = check
        if transfer:
            self.get()

    def get(self, check=2, ifile=None):
        """Download remote files when needed

        :Params: *check*: checks

            - ``0``: check nothing => force the transfer
            - ``1``: only check the existence of the local file
            - ``2``: check existence and compare local and remote dates
        """
        if check is None: check = self.check
        if check is None: check = 2

        # Direct access
        if self.cpmode==0:
            locfiles = self.local_files(ifile)
            if not locfiles:
                self.error('No input file found with this pattern: '+self.local_pattern)
            return locfiles

        # Copy/transfer
        locfiles = []
        tag = 'Downloading %s:'%self.host if self.cpmode==2 else 'Copying '
        for i, remfile in enumerate(self.remote_files(ifile)):

            ## Cache names
            #self.rem2loc[remfile] = locfile
            #self.loc2rem[locfile] = remfile

            # Theoretical local file name
            locfile = os.path.join(self.local_dir, remfile[len(self.remote_dir)+1:])

            # Check
            nothere = not os.path.exists(locfile)
            if check==1 or (check and nothere):
                getit = nothere
                msg = 'no local copy'
            elif check==2:
                loctime = os.stat(locfile)[8]
                if self.cpmode==1:
                    remtime = os.stat(remfile)[8]
                else:
                    remtime = self.ftp.stat(remfile).st_mtime
                getit = remtime > loctime
                msg = 'local copy too old'
            else:
                getit = 1
                msg = 'force transfer'

            # Download
            if getit:

                # Check local dir and file
                locdir = os.path.dirname(locfile)
                if not os.path.exists(locdir):
                    os.makedirs(locdir, mode=self.dmode)
                    self.debug('Create local directory: '+locdir)
                if os.path.exists(locfile):
                    os.remove(locfile)

                # Transfer
                self.debug('%s%s to %s (%s)'%(tag, remfile, locfile, msg))
                if self.cpmode==1: # local
                    shutil.copy2(remfile, locfile)
                else: # remote
                    #self.ftp.get(remfile, locfile)
                    subprocess.check_call(['scp', '%s:%s'%(self.host, remfile), locfile],stdout=subprocess.PIPE)
                os.chmod(locfile, self.fmode)

            locfiles.append(locfile)
        if not locfiles:
            self.error('No file found on remote host with current pattern (%s) and filter'%self.remote_pattern)
        return locfiles

    def remote_files(self, ifile=None):
        """List of remote files"""
        # List
        if self.cpmode<=1: # local copy
            files = glob.glob(self.remote_pattern)
        else: #remote copy
            files = [os.path.join(self.remote_dir, remfile)
                for remfile in self.remote_exec('ls '+self.remote_pattern)]

        # Selection
        if not len(files):
            self.error('No remote files')
        else:
            if callable(self.filter):
                files = self.filter(files)
            if self.sort is not None:
                files.sort(self.sort if self.sort is not True else None)
            if isinstance(ifile, int): return files[ifile]
        return files

    def local_files(self, ifile=None, update=None):
        """List of local (working) files"""
        # List
        if self.cpmode==0:
            files = glob.glob(self.local_pattern)
            if callable(self.filter):
                files = self.filter(files)
            if self.sort is not None:
                files.sort(self.sort if self.sort is not True else None)
        elif self._local_files is None or update is not False: # Update copy
            files = self.get()
            self._local_files = files
        else: # Direct access
            files = self._local_files
        self._local_files = files


        # Selection
        if not len(files):
            self.warning('No local files')
        else:
            if isinstance(ifile, int): return files[ifile]
        return files


class OutputWorkFile(WorkFile):
    """A class to deal with an output remote file

    :Params:

        - *path*: path in the form ``"<prefix>(<remote_dir>><local_dir>)<pattern>"``
           (see :func:`parse_path`)

    :Params:

        - *logger*: a :class:`~vacumm.misc.io.Logger` (or subclass) instance
        - *sshbank*: a :class`SSHBank` instance

    :Example:

    >>> wfile = OutputWorkFile('(sftp://caparmor-sftp/home125>/home200/caparmor)toto/data/file.png')
    >>> wfile = OutputWorkFile('(/home200/caparmor>sftp://username@my.host.fr:1022/prefix)/toto/data/file.png')
    >>> wfile = OutputWorkFile('/home510/toto/toto.png') # does nothing !
    >>> wfile = OutputWorkFile('(/home15>/home12)toto/toto.png') # local copy only
    >>> pylab.savefig(wfile.local_file)
    >>> wfile.put() # send or copy

"""
    def __init__(self, path, rootdir='', logger=None, sshbank=None, subst=None, raise_error=True, **kwargs):
        self._raise_error = raise_error
        # Parse path
        locdir,remdir,basefile = self.parse_path(path, subst=subst, rootdir=rootdir)
        if locdir is None:
            self.host = ''
            self.local_file = self.remote_file = basefile
            self.cpmode = 0
            WorkFile.__init__(self, logger, False, raise_error=raise_error, **kwargs)
            locdir = os.path.dirname(self.local_file)
            if not os.path.exists(locdir):
                os.makedirs(locdir, self.dmode) # local dir
            return
        self.local_file = os.path.join(locdir, basefile)

        # Extract host and port
        up = urlparse(remdir)
        if up.hostname is None:
            self.host = ''
            self.cpmode = 1
        else:
            self.host = up.netloc
            self.cpmode = 2
        remdir = up.path
        self.remote_file = os.path.join(remdir, basefile)

        # Get SSH agent from bak
        ssh = self._host2ssh_(sshbank)

        # Final init
        WorkFile.__init__(self, logger, ssh, raise_error=raise_error, **kwargs)

        # Local dir
        locdir = os.path.dirname(os.path.join(locdir,basefile))
        if not os.path.exists(locdir):
            os.makedirs(locdir, self.dmode)


    def put(self, checkdir=True):
        """Send the file"""
        if self.cpmode == 0:
            return self.local_file
        if not os.path.exists(self.local_file):
            self.warning('Local file to send not found: '+self.local_file)
            return

        # Check remote directory
        if checkdir:
            remdir = os.path.dirname(self.remote_file)
            if self.cpmode==1: # local
                if not os.path.exists(remdir):
                    os.makedirs(remdir, self.dmode)

            else: # remote
                try: # first, check full path
                    self.ftp.stat(remdir)

                except: # now check all subdirs
                    remdirs = remdir.split(os.path.sep)
                    for i in xrange(1, len(remdirs)):
                        thisdir = os.path.join(os.path.sep, *remdirs[:i+1])
                        try:
                            self.ftp.stat(thisdir)
                            continue
                        except:
                            try:
                                self.ftp.mkdir(thisdir)
                                self.ftp.chmod(thisdir, self.dmode)
                            except:
                                self.warning('Cannot create remote dir for transfer: %s:%s'%(self.host, thisdir))
                                return self.local_file

        # Transfer
        if self.cpmode==1: # local
            self.debug('Copying %s to %s'%(self.local_file, self.remote_file))
            try:
                shutil.copy(self.local_file, self.remote_file)
                os.chmod(self.remote_file, self.fmode)
            except:
                self.warning("Can't copy file")
        else: # remote
            remfile = self.remote_file
            if not os.path.isabs(remfile): remfile = '/'+remfile
            self.debug('Uploading %s to sftp://%s%s'%(self.local_file, self.host, remfile))
            try:
                self.ftp.remove(self.remote_file)
            except:
                try:
                    self.ftp.stat(remfile)
                except:
                    pass # File does not exist so no warning
                else:
                    self.warning("Can't remove remote file")
            try:
                self.ftp.put(self.local_file, self.remote_file)
            except:
                self.warning("Can't upload file")
            try:
                self.ftp.chmod(self.remote_file, self.fmode)
            except:
                self.warning("Can't change permission of remote file")
        return self.local_file

