.. highlight:: bash

.. _user.install.installer:

Conda installer script
======================

Here is an example of bash script with options
that can help you install :program:`conda`, vacumm and its dependencies in
a single command.
You can either execute this script or create yours, inspired from this one.

You can download the script :download:`here <../../../scripts/install/install-vacumm-conda-full.sh>`.

Usage, with defaults depending on platform::

    Usage:

    install-vacumm-conda-full.sh [-h|--help] [{-p|--prefix} PREFIX] [{-c|--condainstaller} CONDA_INSTALLER]

     PREFIX: installation prefix (default: $HOME/miniconda2)
     CONDA_INSTALLER: conda installer full path (default: $HOME/Miniconda2-latest-Linux-x86_64.sh)


Content:

.. literalinclude:: ../../../scripts/install/install-vacumm-conda-full.sh
