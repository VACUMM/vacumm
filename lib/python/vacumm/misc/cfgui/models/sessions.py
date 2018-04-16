#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2017)
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

import os

from vacumm.misc.config import ConfigManager, ConfigObj

from vacumm.misc.cfgui import XmlConfigObject, XmlConfigDictObject


class Session(XmlConfigObject):
    '''
    A session managing a configuration and its specification
    '''

    xml_attributes = {
        'name':(None,''),
    }
    xml_textnodes = {
        'specification_file':{'single':True, 'default':''},
        'configuration_file':{'single':True, 'default':''},
    }

    specification_extension = '.ini'
    specification_file_filter = 'Specification files (*%(specification_extension)s) (*%(specification_extension)s);;All files (*) (*)'%locals()

    configuration_extension = '.cfg'
    configuration_file_filter = 'Configuration files (*%(configuration_extension)s) (*%(configuration_extension)s);;All files (*) (*)'%locals()


    def __init__(self, init=None, **kwargs):
        '''
        :Params:
            - **init**: call initilize, default is True
        '''
        XmlConfigObject.__init__(self, **kwargs)
        self._prevspcfile = None
        self._prevcfgfile = None
        self.session_config_manager = ConfigManager()
        self.session_config = ConfigObj()
        # ConfigManager.load apply defaults, session_config_raw is used 
        # to know which options were present or not in the configuration file
        self.session_config_raw = ConfigObj()
        if init: self.initialize()


    session_configspec = property(
        fget=lambda self: self.session_config_manager._configspec if self.session_config_manager else None,
        doc='The internal :class:`ConfigObj` configuration specification')


    def initialize(self, reload=True, reloadspc=None, reloadcfg=None, patch=None):
        '''
        :Params:
            - **reload**: reload everything even if files didn't changed, default is True
            - **reloadspc**: reload specification even if specification file didn't changed, defaults to reload argument value.
            - **reloadcfg**: reload configuration even if configuration file didn't changed, defaults to reload argument value.
        '''
        if reloadspc is None: reloadspc = reload
        if reloadcfg is None: reloadcfg = reload
        self.logger.verbose('Initialize, reloadspc: %(reloadspc)s, reloadcfg: %(reloadcfg)s, patch: %(patch)s', locals())

        if not self.session_config_manager or reloadspc:
            self.logger.verbose('Loading specification: %r', self.specification_file)
            self.session_config_manager = ConfigManager(self.specification_file, interpolation=False)
            self._prevspcfile = self.specification_file
        self.logger.debug('Specification (%r):\n%s\n', self.specification_file, self.pformat(self.session_configspec.dict()))

        if not self.session_config or reloadcfg:
            self.logger.verbose('Loading configuration: %r', self.configuration_file)
            if self.configuration_file:
                if not os.path.exists(self.configuration_file):
                    self.logger.warning('Configuration file does not exists: %r', self.configuration_file)
                self.session_config = self.session_config_manager.load(self.configuration_file)
                self.session_config_raw = ConfigObj(self.configuration_file)
            else:
                self.session_config = self.session_config_manager.defaults()
                self.session_config_raw = ConfigObj()
            self._prevcfgfile = self.configuration_file

        if patch:
            self.logger.debug('Patching with configuration:\n%r\n', self.pformat(patch.dict()))
            self.session_config_manager.patch(self.session_config, patch)

        self.logger.debug('Configuration (%r):\n%s\n', self.configuration_file, self.pformat(self.session_config.dict()))
        return self


    def info(self):
        self.logger.notice('Session:\n%s', self.to_xml_str())
        self.logger.verbose('Specification (%r):\n%s\n', self.specification_file, self.pformat(self.session_configspec.dict()) if self.session_configspec else None)
        self.logger.verbose('Configuration (%r):\n%s\n', self.configuration_file, self.pformat(self.session_config.dict()) if self.session_config else None)


class Sessions(XmlConfigDictObject):

    xml_childnodes = {'sessions':{'type':Session, 'key':'name', 'self':True}}

