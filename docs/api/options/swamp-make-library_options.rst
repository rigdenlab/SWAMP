.. _swamp_make_library_options:

swamp-make-library options
---------------------------

.. code-block:: text

    usage: swamp-make-library    [-h] [-nprocs [NPROCS]] [-homologs [FILENAME]]
                                 [-overwrite_library] [-min_samples [MIN_SAMPLES]]
                                 [-xi [XI]] [-cluster_method [CLUSTER_METHOD]]
                                 [-max_eps [MAX_EPS]] [-eps [EPS]]
                                 [-min_cluster_size [MIN_CLUSTER_SIZE]]
                                 workdir

positional arguments:
+++++++++++++++++++++

.. code-block:: text


      workdir               Working directory for SWAMP-MAKE-LIBRARY.


optional arguments:
+++++++++++++++++++

.. code-block:: text

  -h, --help            show this help message and exit

  -nprocs [NPROCS]      Number of processors to use

  -homologs [FILENAME]  A file with the list of homolog structures to exclude
                        from the resulting library

  -overwrite_library    If set, overwrite the SWAMP library with the new
                        ensembles

  -min_samples [INT]
                        sklearn.OPTICS: no. of samples in a neighborhood for a
                        point to be considered as a core

  -xi [XI]              sklearn.OPTICS: min.steepness on the reachability plot
                        to constitute a cluster boundary

  -cluster_method [CLUSTER_METHOD]
                        sklearn.OPTICS: extraction method using the calculated
                        cluster reachability

  -max_eps [MAX_EPS]
                        sklearn.OPTICS: max. dist. between points to consider
                        within neighborhood of each other

  -eps [EPS]            sklearn.OPTICS: max. dist. between points to consider
                        within neighborhood of each other

  -min_cluster_size [MIN_CLUSTER_SIZE]
                        sklearn.OPTICS: min. no. of samples in a cluster


