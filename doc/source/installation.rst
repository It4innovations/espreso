
===================================
Installation of the library ESPRESO
===================================

Download the library
--------------------

Jak stahnout .. asi stable vetev na git labu


Directory structure
^^^^^^^^^^^^^^^^^^^

zakladni popis adresarove struktury


Dependencies
------------

co vsechno je nutne pro beh ESPRESA


Building the ESPRESO
--------------------

For compilation and instalation ESPRESO uses a Python-based framework `Waf <https://waf.io/book/>`_.
The compilation includes two phases: ``configuration`` and ``installation``.
The former configures persistend data and checks all required headers and libraries for your installation.
The latter makes the library and produces the runable binary file ``espreso``.

Configuration
^^^^^^^^^^^^^

The **configuration** is made by the command: ::

  $ ./waf configure

It configures the library based on the default settings from the file ``build.config.default``.
When something goes wrong, you may want to change some attributes. It should not be done
directly in the file ``build.config.default``. It is recommended to create new file ``build.config``
from the default settings: ::

  $ cp build.config.default build.config

When you create the file ``build.config``, all setttings from this file rewrite attributes
from the default file.

Support for various direct solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MKL, PARDISO, MUMPS


Hardware Acceleration
^^^^^^^^^^^^^^^^^^^^^
CPU, GPU, MIC


Installation
^^^^^^^^^^^^

waf install


Set up the environment
----------------------

cesty + pocet vlaken atd


Testing the installation
------------------------

jak rozbehnout testy
