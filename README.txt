VACUMM
======

VACUMM is primarily designed for the validation of ocean models.
It is made of a python library and a few scripts that help pre-process
and post-process oceanic and atmospheric data from models and observations.
It has some specialized modules for managing model outputs and making
advanced diagnostics.


Requirements
------------

- UVCDAT wich includes other requirements: python 2.7+numpy, 
  scipy, matplotlib, basemap
- configobj

Needed for some modules:

- paramiko
- seawater

See also: http://www.ifremer.fr/vacumm/user.install.prereq.html


Copyright
---------

VACUMM is under the CeCiLL license, which is compatible with well knwon 
GPL license. 

See also: http://www.ifremer.fr/vacumm/appendix.license.html


Documentation
-------------

The documentation explain how to install the librairie, provides numerous 
examples of its usage and describes its full API.

See: 

- official web site: http://www.ifremer.fr/vacumm
- rebuilt everyday: http://relay.actimar.fr/vacumm
- complete: http://relay.actimar.fr/~raynaud/vacumm



Install
-------

1. Install all requirements.
2. Download it following instructions of this page: 
   http://www.ifremer.fr/vacumm/user.install.download.html.
   Then go to the source directory.
3. Run the standard setup script::

    python setup.py install
   
   More info on this page:
   http://www.ifremer.fr/vacumm/user.install.installs.html
   

Contacts
--------

Stephane Raynaud (raynaud (at) actimar.fr)
Guillaume Charria (Guillaume.Charria (at) ifremer.fr)
