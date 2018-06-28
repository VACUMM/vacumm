.. _user.scripts.showmpi:

:program:`showmpi.py` -- Display the MPI decomposition domain
=============================================================

MPI decomposition information comes from a file tipically names :file:`mpi.txt`.

.. include:: bin/showmpi.help.txt

:Examples:

    .. code-block:: bash

        $> showmpi.py -o showmpi.png ../data/mpi.nc ../data/mpi.txt

    .. figure:: images/showmpi.png

        MPI blocks are plotted along indices over the model mask.
        Contour of MPI blocks containing inactive ocean points are in red.
        Inactive ocean points are also colored in red.



:See also: :ref:`user.scripts.showgrid`, :class:`~vacumm.data.mode.mars.Mars3D`.