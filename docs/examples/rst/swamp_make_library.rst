.. _swamp_make_library:

Make a helical pair ensemble library with SWAMP
-----------------------------------------------

It is possible to alter the default SWAMP library using ``swamp-make-library`` a command that allows the user to create a new library of helical pairs and overwrite the original one.


1. Running swamp-make-library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user can create a SWAMP library by simply running the following the command:

.. code-block:: shell

    swamp-make-library "<workdir>"

SWAMP will then save all the ensembles of the library into "<workdir>". The behaviour of the command can be set by using the additional available options, check them out :ref:`here <swamp_make_library_options>`. Most of the optional arguments control how the clusters are created, anf they are directly passed to `scikit-learn OPTICS clustering module <https://scikit-learn.org/stable/modules/generated/sklearn.cluster.OPTICS.html#sklearn.cluster.OPTICS>`_. This command is specially useful if you think SWAMP would do a better job for your target structure with more compact/loose ensembles, or you want to remove homolog structures from the library for benchmarking purposes.
