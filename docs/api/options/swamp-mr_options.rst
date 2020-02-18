.. _swamp_mr_options:

swamp-mr options
----------------

.. code-block:: text

usage: swamp-mr    [-h] [-nprocs [NPROCS]] [-pdb_benchmark [PDB_BENCHMARK]]
                   [-platform [PLATFORM]] [-environment [ENVIRONMENT]]
                   [-queue [QUEUE]] [-mtz_phases [MTZ_PHASES]]
                   [-job_kill_time [JOB_KILL_TIME]]
                   [-ncontacts_threshold [NCONTACTS_THRESHOLD]]
                   [-consco_threshold [CONSCO_THRESHOLD]]
                   [-python_interpreter [PYTHON_INTERPRETER]]
                   [-combine_searchmodels] [-use_centroids] [-use_cores]
                   id workdir mtzfile fastafile conpred sspred



positional arguments:
+++++++++++++++++++++

.. code-block:: text

      id                    Unique identifier for this MR run

      workdir               Working directory for SWAMP-MR

      mtzfile               MTZ file with the reflection data

      fastafile             FASTA file with the sequence of the structure

      conpred               Residue contact prediction for the target structure
                            (psicov format)

      sspred                Secondary structure prediction for the target protein
                            (topcons format)

optional arguments:
+++++++++++++++++++

.. code-block:: text

      -h, --help            show this help message and exit

      -nprocs               Number of parallel processors to use (default 1)

      -pdb_benchmark        PDB file with the solved structure (for benchmarking)

      -platform             Platform to execute MR jobs ('local', 'sge', 'slurm',
                            'local')

      -mtz_phases           MTZ file with phase information (for benchmarking)

      -job_kill_time        Maximum runtime of each MR job

      -ncontacts_threshold  Minimum no. of interhelical contacts to compute CMO of
                            a subtarget against swamp library.

      -consco_threshold     Minimum CMO between predicted and observed contacts to
                            use search model.

      -python_interpreter   Indicate a python interpreter for MR jobs.

      -combine_searchmodels If set, combine search models matching different
                            subtargets

      -use_centroids        Ensemble centroids used as search models as well
                            (only supported if not combining search models)

      -use_cores            Core ensembles used as search models as well
                            (only supported if not combining search models)
