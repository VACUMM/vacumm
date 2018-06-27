.. _user.install.prereqalt:

Prerequisites - alternative installation
========================================

**This page try to describe a "non-official" vacumm's dependencies installation.**

A ligther installation can be achieved using the core required packages.

This kind of installation allows the users to install vacumm with some packages
available in their (recent) OS repositories or by dowloading and installing them
manually ("python setup.py install --user" or "easy_install --user" for example).

Required dependencies
---------------------

Install each package with respect to the minimal (and possibly maximal) required version followed by
the vacumm pacakge itself as describe later in this documentation.

- numpy (1.5.1)
- scipy (0.10.0)
- cdat_lite (6.0rc2) **see notes below**
- matplotlib (1.1.0)
- basemap (1.0.6)
- configobj (4.7.2)


Optionnal dependencies
----------------------

Please refer to the :ref:`user.install.prereq` page for other optionnal packages.

Notes
-----

cdat_lite
^^^^^^^^^

`cdat_lite <http://proj.badc.rl.ac.uk/cedaservices/wiki/CdatLite>`_ is a subset
of CDAT containing the mostly used packages (cdms2, MV2, cdtime, ...)

Vacumm should work with the latest cdat_lite version, we recommend the use of the
latest git version using:

git clone http://proj.badc.rl.ac.uk/git/cdat_lite

**(tested with the revision 557a0b50609e77f75c05b78a15206b392f51dac9 commited the 23/11/2012)**

Some archives are available but not all of vacumm's features may work properly
with the latest (cdat_lite-6.0rc2)


