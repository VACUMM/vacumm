.. _user.tut.misc.file:

File utilities
==============

See also the api documentation: :mod:`~vacumm.misc.file`


Introduction
~~~~~~~~~~~~

The :mod:`~vacumm.misc.file` module provides various file system related features such as:
    - filesystem traversal with depth support
    - file search, wildcard or regex based
    - file rollover (backup)
    - size parsing and formatting
    - directory creation without error on existing directory


File system traversal with depth support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`~vacumm.misc.file.walk`


File search
~~~~~~~~~~~

The find/xfind/efind/xefind functions are intended to provide flexible file or 
directories lookup.

xfind and xefind are generator functions their call do not build a potentially
huge file/directory list so that you can process the result in a loop or use the
optionnal callbacks described below.

find and efind simply call xfind and xefind to return the lokkup result as a list.

The simpliest way to make a search is done by using wildcard pattern (?*[seq][!seq])
with find/xfind. See :mod:`fnmatch` for wildcard patterns details, how they work and
are supported.

Sometimes the need for a more precise lookup appear, the regular expression are a
very powerfull way for that goal although they need a bit of understanding. See
:mod:`re` and http://www.expreg.com/.
Another advantage to the regex method is that you can extract parts of the file path
you selected by using the getmatch parameters with efind/xefind.

find and efind have nearly common features but few differences exists:
    - find/xfind support single or list of pattern/exclude, efind/xefind do not,
      this should not be needed because of the regular expression capabilities.
    - efind/xefind support returning (matching path, match object) couples, find/xfind do not

Using wildcard expressions: See :func:`~vacumm.misc.file.xfind` and :func:`~vacumm.misc.file.find`

Using regular expressions: See :func:`~vacumm.misc.file.xefind` and :func:`~vacumm.misc.file.efind`

The example
-----------

.. literalinclude:: ../../../../scripts/tutorials/misc.file.find.py
    :language: python

.. program-output:: ../../../scripts/tutorials/misc.file.find.py


Making file backups (called rollover)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See :func:`~vacumm.misc.file.rollover`

.. literalinclude:: ../../../../scripts/tutorials/misc.file.rollover.py
    :language: python

.. .program-output:: ../../../scripts/tutorials/misc.file.rollover.py


Displaying and analysing file sizes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See :func:`~vacumm.misc.file.strfsize` and :func:`~vacumm.misc.file.strpsize`

.. literalinclude:: ../../../../scripts/tutorials/misc.file.strfsize.py
    :language: python

.. program-output:: ../../../scripts/tutorials/misc.file.strfsize.py


Creating directories
~~~~~~~~~~~~~~~~~~~~

See :func:`~vacumm.misc.file.mkdirs` and :func:`~vacumm.misc.file.mkfdirs`

.. literalinclude:: ../../../../scripts/tutorials/misc.file.mkdirs.py
    :language: python

.. program-output:: ../../../scripts/tutorials/misc.file.mkdirs.py


