.. _user.faq.misc:

Misc
====


.. _user.faq.misc.batch:

How to submit a script in batch mode?
-------------------------------------

Please see how to use the script :ref:`user.scripts.batch`.


.. _user.faq.misc.remote:

How to work with remote files?
------------------------------

.. warning::

    You must install `paramiko <http://www.lag.net/paramiko/>`_ ::

        shell> easy_install paramiko

Thanks to the :mod:`vacumm.misc.remote` module, you can easly work with
remote files.
These files can:

    - in another directory of the same system,
    - on antoher computer you can connect to with SSH.

In addition, a remote file can be:

    - an ensemble of **input** files, like outputs from model ;
    - a single **output** file, like a figure file.


For instance, in case of a list of input files::

    >>> wf = InputWorkfiles('(caparmor-sftp:/home125>/home200/caparmor)toto*/data/file*.nc')
    >>> print wf.remote_files() # targte files
    >>> print wf.local_files() # local files
    >>> print wf.local_files(update=True) # update the list
    >>> for ncfile in wf.local_files(): ....

In this example, the following files are considered:

    - On caparmor-sftp: :file:`/home125/toto*/data/file*.nc`
    - Locally: :file:`/home200/caparmor/toto*/data/file*.nc`

.. note::

    Remote input files are downloaded (or copied) by default only
    the local copy is not present or older.

In output mode (single file)::

    >>> wf = OutputWorkFile('(/home200/caparmor>sftp://caparmor-sftp/home125)toto/data/file.png')
    >>> pylab.savefig(wf.local_file) # write local file
    >>> wf.put() # send it to the remote host

For more information, please the documentation of
the :mod:`vacumm.misc.remote` module.
