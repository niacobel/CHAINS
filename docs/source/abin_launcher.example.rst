**********
Sample run
**********

This section will cover a sample run of ``ABIN LAUNCHER``, to get a concrete example of its operating mode. Our study case here includes three molecules (or geometry files) and two configuration files. Our goal is to optimize the geometry of those three molecules with two different basis sets through ORCA_. 

.. Tip::

   Every file presented in this section can be downloaded `here <https://github.com/niacobel/CHAINS/tree/master/docs/source/sample_files>`_. The other files can be fetched directly from ``ABIN LAUNCHER``'s :ref:`directory <abin_directory>`. Note that for this example, the ``renderer.py`` file from the ``sample_files`` directory replaces the one present in ``ABIN LAUNCHER``'s directory (which is a more complex one, used by CHAINS).

Preparation
===========

Jinja templates
---------------

We will simply use the Jinja templates defined in the :doc:`rendering section <abin_launcher.rendering>` and place them in the ``templates`` directory of ``ABIN LAUNCHER``:

- The Jinja template for the ORCA input file:

.. literalinclude:: sample_files/sample_orca.inp.jinja
   :language: jinja
   :caption: sample_orca.inp.jinja

- The Jinja template for the job script:

.. literalinclude:: sample_files/sample_orca_job.sh.jinja
   :language: jinja
   :caption: sample_orca_job.sh.jinja

Rendering function
------------------

We will also use the rendering function we defined in the :ref:`rendering section <rendering_fct>`, this function is placed in the ``renderer.py`` file. For our sample case, this file contains:

.. literalinclude:: sample_files/renderer.py
   :language: python
   :caption: renderer.py

Clusters configuration file
---------------------------

For this run, we will be running on the LEMAITRE3 cluster from the CECI_, the job scheduler is SLURM and the chosen scaling function is ``total_nb_elec``. 

The clusters configuration file looks like:

.. literalinclude:: sample_files/clusters.yml
   :language: yaml
   :caption: clusters.yml

Execution
=========

Geometries and configuration files
----------------------------------

Our three molecules are simply |CH4|, |C2H6| and |C3H8|. They are each represented by a XYZ geometry file:

.. |CH4| replace:: CH\ :sub:`4`\ 
.. |C2H6| replace:: C\ :sub:`2`\ H\ :sub:`6`\ 
.. |C3H8| replace:: C\ :sub:`3`\ H\ :sub:`8`\ 

.. literalinclude:: sample_files/ch4.xyz
   :caption: ch4.xyz

.. literalinclude:: sample_files/c2h6.xyz
   :caption: c2h6.xyz

.. literalinclude:: sample_files/c3h8.xyz
   :caption: c3h8.xyz

For our two configuration files, we want to perform a geometry optimization using the DFT method with the B3LYP functional, and the two considered basis sets are def2-SVP and def2-TZVP. Thus, our two files are:

.. literalinclude:: sample_files/svp.yml
   :language: yaml
   :caption: svp.yml

.. literalinclude:: sample_files/tzvp.yml
   :language: yaml
   :caption: tzvp.yml

Let's store those files in two distinct directories: one for the molecule geometries, named ``molecules``, and one for the configuration files, named ``configs``. We will also create a directory for the jobs, named ``orca_jobs``.

This is what our directory on the cluster, named ``abin_docs_sample``, now looks like:

.. code-block::

   abin_docs_sample/
      └── abin_launcher/ 
            ├── abin_launcher.py
            ├── geom_scan.py
            ├── scaling_fcts.py
            ├── renderer.py
            ├── abin_errors.py
            ├── clusters.yml
            ├── mendeleev.yml
            └── templates/
                  ├── sample_orca.inp.jinja
                  └── sample_orca_job.sh.jinja
      └── molecules/ 
            ├── ch4.xyz
            ├── c2h6.xyz
            └── c3h8.xyz
      └── configs/ 
            ├── svp.yml
            └── tzvp.yml
      └── orca_jobs/ 
            └── currently empty

Running ABIN LAUNCHER
---------------------

.. Tip::

   Before executing ``ABIN LAUNCHER``, remember to load (manually or through your user profile configuration) your Python distribution (version 3.5+), which must include PyYAML (version 5.1+) and Jinja2 (version 2.10+).

We can now execute ``abin_launcher.py`` by running the command (from ``abin_docs_sample``):

.. code-block:: console

   $ python abin_launcher/abin_launcher.py -m molecules/ -cf configs/ -p orca -o orca_jobs/ -cl lemaitre3

This is what appears on the console screen (here captured in a file for ease of reading):

.. literalinclude:: sample_files/sample_run.log
   :language: text
   :caption: sample_run.log

And if we take a look at the job queue, we can see that our 6 ORCA jobs have indeed been submitted to the job scheduler:

.. figure:: figures/abin_launcher_squeue.*
    :scale: 70%
    :align: center
    :alt: Jobs queue after running ABIN LAUNCHER
    :figclass: align-center
    
    Jobs queue after running ``ABIN LAUNCHER``

Directory structure after execution
-----------------------------------

If we take a look at the directory structure, it now looks like:

.. code-block::

   abin_docs_sample/
      └── abin_launcher/ 
            └── no changes
      └── molecules/
            └── launched/
                  ├── ch4.xyz
                  ├── c2h6.xyz
                  └── c3h8.xyz
      └── configs/ 
            └── launched/
                  ├── svp.yml
                  └── tzvp.yml
      └── orca_jobs/ 
            └── ch4_svp/
                  ├── ch4.xyz
                  ├── svp.yml
                  ├── ch4.inp
                  ├── orca_job.sh
                  └── ch4_svp.log
            └── ch4_tzvp/
                  ├── ch4.xyz
                  ├── tzvp.yml
                  ├── ch4.inp
                  ├── orca_job.sh
                  └── ch4_tzvp.log
            └── c2h6_svp/
                  ├── c2h6.xyz
                  ├── svp.yml
                  ├── c2h6.inp
                  ├── orca_job.sh
                  └── c2h6_svp.log
            └── c2h6_tzvp/
                  ├── c2h6.xyz
                  ├── tzvp.yml
                  ├── c2h6.inp
                  ├── orca_job.sh
                  └── c2h6_tzvp.log
            └── c3h8_svp/
                  ├── c3h8.xyz
                  ├── svp.yml
                  ├── c3h8.inp
                  ├── orca_job.sh
                  └── c3h8_svp.log
            └── c3h8_tzvp/
                  ├── c3h8.xyz
                  ├── tzvp.yml
                  ├── c3h8.inp
                  ├── orca_job.sh
                  └── c3h8_tzvp.log

As you can see, since they have been successfully processed, the geometry files and the configuration files have both been archived into a ``launched`` directory created by ``ABIN LAUNCHER``. This is the default behavior, allowing you to reuse the same directories for other geometries and configurations, making it easier for example to create an alias for the execution command. If you want to turn it off, just add the ``-km / --keep_mol`` and/or ``-kc / --keep_cf`` optional arguments to keep the geometry files and/or the configuration files, respectively.

A subdirectory has also been created in ``orca_jobs`` for each of the six jobs. Those subdirectories contain the copies of the geometry and the configuration files, the rendered input file and job script, as well as a log file containing the details of the treatment of this geometry-configuration combination by ``ABIN LAUNCHER``. *(The output files created by ORCA will also end up in those subdirectories.)*

As an example, here is what the ``c3h8_tzvp.log`` file looks like:

.. literalinclude:: sample_files/c3h8_tzvp.log
   :language: text
   :caption: c3h8_tzvp.log

.. Hyperlink targets

.. _CECI: http://www.ceci-hpc.be/
.. _ORCA: https://www.faccts.de/orca/