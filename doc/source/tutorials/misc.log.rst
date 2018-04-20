.. _user.tut.misc.log:

Logging utilities
=================

See also the api documentation: :mod:`~vacumm.misc.log`

The :class:`vacumm.misc.log.Logger` class provides shortcuts for logging features:
    - Some non standard levels are added (verbose and notice)
    - Logging can be customized with different configurations for messages sent to 
      the standard output (terminal) and to a file.
    - Some configuration can be massivelly changed
    - Handlers can be shared between several loggers

Sample code using loggers:

.. literalinclude:: ../../../../scripts/tutorials/misc.log.py
    :language: python

The code above produce the following output:

.. command-output:: ../../../scripts/tutorials/misc.log.py

