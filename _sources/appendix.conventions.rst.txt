.. _appendix.conventions:

Coding rules
************

Coding norm : PEP 8
===================

The coding rules mainly follow he standards specified by the :pep:`8`.
These standards include the following important points:

    - The indentations are composed of four spaces (no tabs).
    - The names of classes and exceptions are capitalized at the beginning of the word,
      and no character is used to separate words in the case of compound nouns. Example : ``MySuperClass``.
    - Other names (functions, methods, attributes) are all lowercase and
      underscore `` _`` character is used to separate words. Example: ``my_super_function``.
    - The encoding must be UTF-8.

The following points are of secondary importance:

    - Lines must not exceed 79 characters.
    - Spaces must be used to separate operators.


Docstrings
==========

The docstrings are  comments located immediately after a statement.
They are used to learn about the stated purpose (purpose, parameters, outputs ...)
and this information will be inserted in the documentation of the library after compilation.

Docstrings must be written in `rst <http://docutils.sourceforge.net/rst.html>`_ programming language,
used by the documentation generator `Sphinx <http://sphinx.pocoo.org>`_ .
The language used is English.

 We propose the following form::

    def myfunc(a, b, c=None, d=5, **kwargs):
        """Short description of myfunc

        More information about myfunc.
        You can even include images.

        .. note:: A simple note.

        :Params:

            - **a**: First mandatory parameter.
              Its description continues here.
            - **b**: The second mandatory one.
            - **c**, optional: First optional keyword parameter.
            - **d**, optional: Second one.
            - Other parameters.

        :Return:

            - *x*: First return variable.
            - *y*: Second one.

        :Example:

            >>> x, y = myfunc(a, b, d=4, g='name')
        """
        ...

