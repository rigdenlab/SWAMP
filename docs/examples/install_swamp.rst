Installing SWAMP
----------------

1. Ensure dependencies
^^^^^^^^^^^^^^^^^^^^^^^
Before installing SWAMP you may wish to make sure that the following `dependencies <https://github.com/rigdenlab/SWAMP/tree/master/docs/requirements.txt>`_ are installed. Additionally, SWAMP requires a locally installed version of CCP4.

2. Clone the repository
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    git clone https://github.com/rigdenlab/SWAMP

3. Create symbolic links to CCP4 libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    ln -s SWAMP/swamp $CCP4/lib/py2/swamp
    ln -s SWAMP/bin/swamp-mr $CCP4/bin/swamp-mr
    ln -s SWAMP/bin/swamp-results $CCP4/bin/swamp-results
    ln -s SWAMP/bin/swamp-make-library $CCP4/bin/swamp-make-library
