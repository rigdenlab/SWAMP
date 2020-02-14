.. _docs_install:

Installing SWAMP
----------------

1. Ensure dependencies
^^^^^^^^^^^^^^^^^^^^^^^
Before installing SWAMP you may wish to make sure that the following `dependencies <https://github.com/rigdenlab/SWAMP/tree/master/docs/requirements.txt>`_ are installed. Additionally, SWAMP requires a locally installed version of CCP4.

2. Clone the repository
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    git clone https://github.com/rigdenlab/SWAMP

3. Build SWAMP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir SWAMP/build
    cd SWAMP/build
    cmake ..
    make install

