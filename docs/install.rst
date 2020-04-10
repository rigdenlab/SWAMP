.. _docs_install:

Installing SWAMP
----------------

Before installing SWAMP, make sure that your system has a CCP4 7.1 installation. You can download and install CCP4 suite `here <http://www.ccp4.ac.uk/index.php>`_.

1. Clone the repository
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    git clone https://github.com/rigdenlab/SWAMP

2. Install dependencies
^^^^^^^^^^^^^^^^^^^^^^^

You can find a complete list of SWAMP dependencies `here <https://raw.githubusercontent.com/rigdenlab/SWAMP/master/requirements.txt>`_. The easiest way to set up all of them is by simply running:

.. code-block:: shell

    ccp4-python SWAMP/setup.py

3. Build SWAMP
^^^^^^^^^^^^^^

Once all the requirements are met the only thing left ot do before using SWAMP is to build it.

.. code-block:: shell

    mkdir SWAMP/build
    cd SWAMP/build
    cmake ..
    make install

