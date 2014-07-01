.. _user.install.download:

Download
========

You can download VACUMM on the GForge site of IFREMER :
https://forge.ifremer.fr/projects/vacumm


Tarballs
--------

Tarballs are available here : https://forge.ifremer.fr/frs/?group_id=93

.. highlight:: bash

To uncompress a tarball, just use for instance::
    
    $ tar xzf vacumm-2.0-svn1231.tar.gz


Subversion
----------

The development version of the library is available on its subversion repository.

For normal or anonymous users::
    
    $ svn checkout --username anonsvn https://forge.ifremer.fr/svn/vacumm/trunk vacumm
    
with password: ``anonsvn``.
 
 
For developpers::

    $ svn checkout --username username https://forge.ifremer.fr/svn/vacumm/trunk vacumm

To update your version, go at the root and type::

    $ svn update

To download a tag, use 