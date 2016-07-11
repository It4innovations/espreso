

Examples
========

Install the library
-------------------

This section serves the common examples,
how to install library to a particular machine.

CPU version with 32-bit length integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is the default settings for ESPRESO.
Hence, no changes of build scripts are needed. ::

  $ git clone git@code.it4i.cz:mec059/espreso.git
  $ cd espreso
  $ ./waf configure
  $ ./waf install

CPU version with 64-bit length integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ESPRESO is 32-bit solver by default.
To allow 64-bit integer, ``build.config`` has to be edited (See `configuration <installation.html#configuration>`__). ::

  $ git clone git@code.it4i.cz:mec059/espreso.git
  $ cd espreso
  $ cp build.config.default build.config

  change value in build.config to: INT_WIDTH = 64

  $ ./waf configure
  $ ./waf install


Accelerated version on a cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Installing ESPRESO on a cluster with accelerator support includes three steps:
  1. Set appropriate modules
  2. Change build.config script
  3. Install

The following example shows the installation with ``Intex Xeon PHI`` support on ``Salomon`` ::

  connect to an accelerated node
  $ . env/modules.salomon

  $ git clone git@code.it4i.cz:mec059/espreso.git
  $ cd espreso
  $ cp build.config.default build.config

  change value in build.config to: SOLVER = MIC

  $ ./waf configure
  $ ./waf install

.. warning ::
  ``env/modules.salomon`` can be out of date.
  It is recommended to change the modules loaded by this script to the newest version,
  if something goes wrong.


Installing on Cray
^^^^^^^^^^^^^^^^^^

Similarly as on a cluster.
The first step is loading of approriate modules.
The following example shows the installation on ``Sisu`` ::

  $ . env/modules.sisu

  $ git clone git@code.it4i.cz:mec059/espreso.git
  $ cd espreso
  $ ./waf configure --cray
  $ ./waf install


Run the solver
--------------

All examples assume that the environment is set (See `environment settings <installation.html#set-up-the-environment>`__).


