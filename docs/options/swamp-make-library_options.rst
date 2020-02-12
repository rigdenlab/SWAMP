SWAMP-MAKE-LIBRARY options
---------------------------

.. code-block::

    usage: swamp-make-library    [-h] [-nprocs [NPROCS]] [-homologs [FILENAME]]
                                 [-overwrite_library] [-core_ensemble]
                                 workdir


positional arguments:
+++++++++++++++++++++

.. code-block::


      workdir               Working directory for SWAMP-MAKE-LIBRARY.


optional arguments:
+++++++++++++++++++

.. code-block::

  -h, --help            show this help message and exit

  -nprocs [NPROCS]      Number of processors to use

  -homologs [FILENAME]  A file with the list of homolog structures to exclude
                        from the resulting library

  -overwrite_library    If set, overwrite the SWAMP library with the new
                        ensembles

  -core_ensemble        If set, ensembles will be trimmed to the core
                        alignment between the constituent models



