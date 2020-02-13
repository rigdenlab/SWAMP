Performing molecular replacement with SWAMP
-------------------------------------------

1. Running swamp-mr
^^^^^^^^^^^^^^^^^^^

In order to run swamp-mr you will first need to obtain `transmembrane topology <http://topcons.cbr.su.se/>`_ and `residue contact <http://raptorx.uchicago.edu/ContactMap/>`_ predictions for your structure of interest. Once this is done, swamp-mr can be executed as follows:

.. code-block:: shell

    swamp-mr "<id>" "<workdir>" "<mtzfile>" "<fastafile>" "<conpred>" "<sspred>"

Additional options are available, check them out `here <https://github.com/rigdenlab/SWAMP/tree/master/docs/options/swamp-mr_options.rst>`_. After running the above command, SWAMP will create the directory ``<workdir>/swamp-mr``. Within this directory, both a ``swamp_scan`` and a ``swamp_mr`` directories will be created, containing the results obtained for the library CMO scan and the molecular replacement, respectively.


A full swamp-mr run can result quite time consuming, if you wish to check out the results swamp has obtained so far you can use ``swamp-results`` as described `here <https://github.com/rigdenlab/SWAMP/blob/master/docs/examples/swamp-results.rst>`_.

2. Output
^^^^^^^^^

Once SWAMP has finished execution of all MR jobs, a table summarising the obtained results will be printed:

+-------------------+---------+-----+-------------+-------------+------------+------------+-------------+-------------+---------+----------+-------------+----------+
|     SEARCH ID     |   LLG   | TFZ | PHSR_CC_LOC | PHSR_CC_ALL | RFMC_RFREE | RFMC_RFACT | RFMC_CC_LOC | RFMC_CC_ALL | SHXE_CC | SHXE_ACL | IS_EXTENDED | SOLUTION |
+===================+=========+=====+=============+=============+============+============+=============+=============+=========+==========+=============+==========+
|  search_277_run_1 |  41.049 | 8.7 |    0.612    |    0.414    |   0.5745   |   0.5355   |     0.65    |    0.434    |  32.65  |   21.0   |     YES     |    YES   |
+-------------------+---------+-----+-------------+-------------+------------+------------+-------------+-------------+---------+----------+-------------+----------+
|  search_951_run_1 |  82.286 | 6.0 |    0.529    |     0.21    |   0.6018   |   0.5927   |    0.608    |    0.248    |  35.33  |   25.0   |     YES     |    YES   |
+-------------------+---------+-----+-------------+-------------+------------+------------+-------------+-------------+---------+----------+-------------+----------+
|  search_169_run_1 |  71.073 | 5.0 |    0.125    |    0.023    |   0.6215   |   0.5897   |    0.112    |    0.026    |  23.08  |   9.0    |     YES     |    NO    |
+-------------------+---------+-----+-------------+-------------+------------+------------+-------------+-------------+---------+----------+-------------+----------+
| search_1008_run_1 |  90.720 | 6.4 |    0.218    |    0.123    |   0.5793   |   0.5922   |    0.228    |    0.125    |  23.03  |   10.0   |     YES     |    NO    |
+-------------------+---------+-----+-------------+-------------+------------+------------+-------------+-------------+---------+----------+-------------+----------+

You can see how to interpret the results displayed in this table `here <https://github.com/rigdenlab/SWAMP/blob/master/docs/examples/swamp-results.rst>`_.