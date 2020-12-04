********
Overview
********

What is ABIN LAUNCHER?
======================

The **Ab Initio Input Builder and Job Launcher**, named ``ABIN LAUNCHER`` for short, is a script that creates input files for ab initio programs, and then submits (or "launches") the corresponding calculations on a job scheduler. It can operate with one or more molecules (through their respective geometry files) and one or more configurations (the parameters for the ab initio program). ``ABIN LAUNCHER`` is the first main script of CHAINS and is executed twice through it:

- The first time, it is used to build the input files associated with the ORCA_ program, which performs the geometry optimization of the molecule. 
- The second time, it is used to build the input files associated with the the Q-CHEM_ program, which performs the calculation of the different properties of our molecule. 

In both cases, after having created the files, ``ABIN LAUNCHER`` then launches the corresponding job on the cluster.

However, ``ABIN LAUNCHER`` is **completely autonomous** and does not depend on CHAINS. It does not need any files outside the ones present in its own directory. It can be extracted and used to launch independent calculations, and can be very easily adapted to deal with other ab initio programs. As such, this specific part of the documentation only explains how ``ABIN LAUNCHER`` works in and of itself, and its integration into CHAINS is explained :doc:`elsewhere <chains.abin_integration>`.

.. _abin_directory:

Directory structure
===================

The ``ABIN LAUNCHER`` `directory <https://github.com/niacobel/CHAINS/tree/master/abin_launcher>`_ has and must keep the following structure:

.. code-block:: text

    abin_launcher/
      ├── abin_launcher.py
      ├── geom_scan.py
      ├── scaling_fcts.py
      ├── renderer.py
      ├── abin_errors.py
      ├── clusters.yml
      ├── mendeleev.yml
      └── templates/

As for what each file does, everything will be explained in more details in the other sections of this documentation, but here is a short summary:

- ``abin_launcher.py`` is the main script itself, the one that needs to be executed.
- ``geom_scan.py`` is the library of functions that define how to scan the geometry files, i.e. how to read and interpret them.
- ``scaling_fcts`` is the library of functions that define how to determine the job size, by calculating what we defined as the "scale index".
- ``renderer.py`` is the library of functions that define how to render the Jinja templates, i.e. how to create the input files and the job script.
- ``abin_errors.py`` contains all the classes and functions defining :ref:`how to handle errors <abin_errors>`.
- ``clusters.yml`` is the YAML file containing all the information specific to the different clusters, called the :ref:`clusters configuration file <clusters_file>`.
- ``mendeleev.yml`` is a YAML version of Mendeleev's periodic table procured by `AlexGustafsson's molecular-data Github repository`_.
- ``templates`` is the directory containing all the Jinja templates that will be used by ``ABIN LAUNCHER``. 

In CHAINS' case, the ``templates`` directory structure is:

.. code-block:: text

   └── templates/
         ├── orca_job.sh.jinja
         ├── orca.inp.jinja
         ├── qchem_job.sh.jinja
         └── qchem.in.jinja

where

- ``orca.inp.jinja`` and ``qchem.in.jinja`` are the Jinja templates corresponding to the input file for ORCA and Q-CHEM, respectively.
- ``orca_job.sh.jinja`` and ``qchem_job.sh.jinja`` are the Jinja templates corresponding to the job script for ORCA and Q-CHEM, respectively. They contain all the commands that will be executed through the job scheduler to perform the calculation.

.. note::

   There are three supplementary files: ``benchmark.py``, ``cron_benchmark.sh`` and, in the ``templates`` directory, ``benchmark.jinja``. Those files are not part of ``ABIN LAUNCHER`` itself but they introduce a :doc:`benchmarking tool <abin_launcher.benchmark>`.

.. _abin_arguments:

Command line arguments
======================

.. argparse::
   :module: abin_launcher
   :func: parser
   :prog: abin_launcher.py
   :nodescription:

How does it work?
=================

The executable part of ``ABIN LAUNCHER`` is the main script, ``abin_launcher.py``. This is the one that must be called in the command line. The overall procedure follows three main steps: **scanning**, **scaling** and **rendering**, followed by the small **submitting** step. Each of the three main steps will be more thoroughly explained in a dedicated section of this documentation. As such, this subsection will only focus on the global procedure.

.. Important::

   An important file that will be often referenced throughout this documentation is the **YAML clusters configuration file** (``clusters.yml``). Rather than presenting it in its entirety at the beginning, the relevant bits of information will be introduced in the different sections, but you can have a full overview of that file in its :ref:`specific documentation <clusters_file>`.

Input files
-----------

There are two main input files for ``ABIN LAUNCHER``:

- :guilabel:`-m / \\--mol_inp`, the **geometry files**.

   The files presenting the nature and the structure of your molecules. They contain the type and number of the constituting atoms and their respective coordinates.

- :guilabel:`-cf / \\--config`, the **configuration files**.

   The YAML files containing the parameters specific to your calculations and your programs (job type, basis set, etc.). Those files must have the .yml or .yaml extension.

In both cases, you can either indicate a specific file in the command line, or point towards a directory where there are multiple of those files. If you specify multiple input files, ``ABIN LAUNCHER`` will process each geometry-configuration combination. For example, if you have 5 geometry files and 3 configuration files, you will end up with 15 launched jobs on your cluster.

By default, every input file that has been *successfully* processed by ``ABIN LAUNCHER`` will be **archived** in a ``launched`` directory created in the same directory as the input files. This has been designed this way so that you can repeatedly use the same directory as "source" for those input files without repeating jobs. If you want to turn off this behavior, you can use the :guilabel:`-km / \\--keep_mol` and/or :guilabel:`-kc / \\--keep_cf` optional arguments to keep the geometry files and/or the configuration files, respectively. Note that if a problem occurs when processing a particular geometry-configuration combination, both the geometry file and the configuration file will be kept in place.

Other arguments
---------------

There are three other required arguments for executing ``ABIN LAUNCHER``:

- :guilabel:`-p / \\--profile`, the **name of the profile** you want to run jobs with.

   The profile is a label used by ``ABIN LAUNCHER`` to know which information to get from its different files. It defines what type of scaling and rendering will be performed, as well as which program will be run and how. As will be explained further in this documentation, this value must be the same as the one given in the :ref:`clusters configuration file <clusters_file>` and in the name of the :ref:`rendering function <rendering_fct>` *(This is case-sensitive!)*.

- :guilabel:`-cl / \\--cluster_name`, the **name of the cluster** you are running on.

   This value must be the same as the one given in the :ref:`clusters configuration file <clusters_file>`, so that ``ABIN LAUNCHER`` knows what you are referring to. *(This is case-sensitive!)*

.. Tip::

   This argument does not need to be the same name as the actual name of your machine. It is just a label used by ``ABIN LAUNCHER`` to know which information to get from its clusters configuration file.

- :guilabel:`-o / \\--out_dir`, the **output directory** 

   This is the directory where each job subdirectory will be created. Those subdirectories are the ones where the files will be created and from which the jobs will be submitted to the job scheduler.

There are also a number of optional arguments that can be used to adapt to each specific situation. Their description in the :ref:`command line arguments <abin_arguments>` subsection should be self-explanatory.

First step: Scanning
--------------------

``ABIN LAUNCHER`` begins by scanning the geometry file, looking for the chemical formula and the atomic coordinates of the molecule. 

For more details on how this scan is performed, consult the :doc:`abin_launcher.scan` specific documentation.

.. Caution::
   At this time, only the XYZ format is supported for geometry files. However, new formats can be added if the need arises.

Second step: Scaling
--------------------

Based on the information received from the geometry file, ``ABIN LAUNCHER`` attributes a value, called the scale index, to the molecule. This value is then used to evaluate the job scale for that molecule and specify the calculation requirements accordingly (walltime, number of CPUs, memory, etc.). 

For more details on how this scaling process is performed, consult the :doc:`abin_launcher.job_scale` specific documentation.

Third step: Rendering
---------------------

Finally, based on user-defined Jinja templates, ``ABIN LAUNCHER`` creates the input files and the job script associated with our calculation. The content of those files is based on the information from the geometry file and the configuration file. 

For more details on how this whole rendering process is performed, consult the :doc:`abin_launcher.rendering` specific documentation.

.. _submitting_step:

The end step: Submitting
------------------------

Now that everything has been prepared for the job, ``ABIN LAUNCHER`` submits it to the job scheduler. The exact command that will be executed is:

.. code-block:: console

    $ <submit_command> <delay_command> <job script>

where

- ``<submit_command>`` is the command which submits jobs to your job scheduler. In SLURM's case, it is the ``sbatch`` command. This must be indicated in the :ref:`clusters configuration file <clusters_file>`: 

   .. code-block:: yaml

      mycluster:
        submit_command: <submit_command>

   where ``mycluster`` is the name of your cluster (the same that was given as the :guilabel:`-cl / \\--cluster_name` command line argument).

- ``<delay_command>`` is an optional command that can delay the submission of a particular job, which can prove useful if you want to prioritize certain job sizes (consult the :doc:`abin_launcher.job_scale` specific documentation for details). In SLURM's case, this is covered by the ``--begin`` argument.
- ``<job script>`` is the name of the file that will be created through the :doc:`rendering process <abin_launcher.rendering>`. It contains the commands needed by the job scheduler to run the calculation on the cluster.

For example, if we want to run an ORCA calculation on a SLURM cluster, but delay the submission of this job by 60 seconds, the command executed by ``ABIN LAUNCHER`` might look like:

.. code-block:: console

    $ sbatch --begin=now+60 orca_job.sh

Once the job has been submitted, ``ABIN LAUNCHER`` will proceed to the next configuration file with the same geometry. Once all the configuration files have been treated, it will proceed to the next geometry and treat again all the configuration files for that geometry. At the end of the execution, barring any problems, a job will have been launched for each geometry-configuration combination.

.. _out_dir_struct:

Output directory structure
--------------------------

If we have for example 2 geometry files and 2 configuration files, once the execution of ``ABIN LAUNCHER`` has ended, the structure of the output directory (given as the :guilabel:`-o / \\--out_dir` command line argument) might look like:

.. code-block:: text

    out_dir/ 
      └── geometry1_config1/
            ├── geometry1.xyz
            ├── config1.yml
            ├── geometry1_config1.log
            ├── job_script.sh
            └── input_file
      └── geometry1_config2/
            ├── geometry1.xyz
            ├── config2.yml
            ├── geometry1_config2.log
            ├── job_script.sh
            └── input_file
      └── geometry2_config1/
            ├── geometry2.xyz
            ├── config1.yml
            ├── geometry2_config1.log
            ├── job_script.sh
            └── input_file
      └── geometry2_config2/
            ├── geometry2.xyz
            ├── config2.yml
            ├── geometry2_config2.log
            ├── job_script.sh
            └── input_file

where 

- ``geometryX_configX`` is the job subdirectory created by ``ABIN LAUNCHER``, and from which the job will be submitted to the job scheduler.
- ``geometryX.xyz`` and ``configX.yml`` are copies of the geometry file and the configuration file, respectively.
- ``job_script.sh`` and ``input_file`` are the files created by the :doc:`rendering process <abin_launcher.rendering>`.
- ``geometryX_configX.log`` is an output file containing the details of the treatment of this geometry-configuration combination by ``ABIN LAUNCHER`` (the computed scale index, the used job scale, etc.)

.. Hyperlink targets

.. _`AlexGustafsson's molecular-data Github repository`: https://github.com/AlexGustafsson/molecular-data
.. _ORCA: https://www.faccts.de/orca/
.. _Q-CHEM: https://www.q-chem.com/
