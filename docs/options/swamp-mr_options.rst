SWAMP-MR options
----------------

SWAMP-MR: Solving structures With Alpha Membrane Pairs - Molecular Replacement

.. code-block::

    usage: swamp-mr.py [-h] [-nprocs [NPROCS]] [-pdb_benchmark [FILENAME]]
                       [-platform [PLATFORM]] [-mtz_phases [FILENAME]]
                       [-job_kill_time [JOB_KILL_TIME]]
                       [-ncontacts_threshold [NCONTACTS_THRESHOLD]]
                       [-consco_threshold [CONSCO_THRESHOLD]]
                       [-python_interpreter [PYTHON_INTERPRETER]]
                       [-combine_searchmodels] [-use_centroids]
                       id workdir mtzfile fastafile conpred sspred


positional arguments:
+++++++++++++++++++++

.. code-block::

      id                    Unique identifier for this MR subroutine

      workdir               Working directory to perform the MR

      mtzfile               MTZ file with the reflection data

      fastafile             FASTA file with the sequence of the structure

      conpred               Residue contact prediction for the target structure

      sspred                Secondary structure prediction for the target protein

optional arguments:
+++++++++++++++++++

.. code-block::

      -h, --help            show this help message and exit

      -nprocs               Number of parallel processors to use

      -pdb_benchmark        PDB file with the solve structure (for benchmarking)

      -platform             Platform to execute MR runs

      -mtz_phases           MTZ file with phase information (for benchmarking)

      -job_kill_time        Maximum runtime of each MR run

      -ncontacts_threshold  Minimum no. of interhelical contacts to compute CMO of
                            subtarget

      -consco_threshold     Minimum CMO between predicted and observed contacts to
                            use search model

      -python_interpreter   Indicate a python interpreter for MR runs

      -combine_searchmodels If set, combine search models matching different parts
                            of the structure

      -use_centroids        Centroids used as search models as well (only
                            supported if not combining search models)


