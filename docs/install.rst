.. _docs_install:

Installing SWAMP
----------------

Before installing SWAMP, make sure that you systems meets the following dependencies:

    * `CCP4 <http://www.ccp4.ac.uk/index.php>`_ >= 7.1
    * `Biopython <https://github.com/biopython/biopython>`_ >= 1.74
    * `PyJob <https://github.com/fsimkovic/pyjob>`_ >= 0.4
    * `ConKit <https://github.com/rigdenlab/conkit>`_ >= 0.11



1. Clone the repository
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    git clone https://github.com/rigdenlab/SWAMP



2. Build SWAMP
^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir SWAMP/build
    cd SWAMP/build
    cmake ..
    make install

