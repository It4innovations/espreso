

Complete list of ESPRESO input parameters
=========================================

The list contains all available input parameters.
Parameters are displayed with the namespace,
however, ESPRESO recongnizes only a parameter name.
To set a particular parameter use the configuration file.
Parameters in the file have the following form: ::

  PARAMETER = VALUE

Parameters can be also set by command line arguments: ::

  $ ./espreso --PARAMETER=VALUE

It overrides the value set in the configuration file.
Detailed description can be found `here <run.html#settings-input-parameters>`_.

Enum types contains allowed values for an appropriate parameter.
For parameter ``PARAM``, the corresponding enum has name ``PARAMalternatives``.
To set a parameter to a particular alternative,
use a number of a particular enum name.


Mesh configuration
------------------

.. doxygennamespace:: espreso::config::mesh
   :members:
   :content-only:

Assembler configuration
-----------------------

.. doxygennamespace:: espreso::config::assembler
   :members:
   :content-only:

Solver configuration
--------------------

.. doxygennamespace:: espreso::config::solver
   :members:
   :content-only:

Output configuration
--------------------

.. doxygennamespace:: espreso::config::output
   :members:
   :content-only:

Debugging configuration
-----------------------

.. doxygennamespace:: espreso::config::info
   :members:
   :content-only:


