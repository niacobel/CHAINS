********
Overview
********

What is CONTROL LAUNCHER?
=========================

The **QOCT-RA Input Builder and Job Launcher**, named ``CONTROL LAUNCHER``, is a script that creates input files for the QOCT-RA_ scripts, and then submits (or "launches") the corresponding calculations on a job scheduler. For a given source file (that describes the needed molecular properties), it can operate with one or more set of QOCT-RA parameters (the configuration files). ``CONTROL LAUNCHER`` is the second main script of CHAINS and is used to execute the second phase of CHAINS' methodology, which consists in evaluating the efficiency of the control procedure and the controllability of the electrons in the molecule.  

However, ``CONTROL LAUNCHER`` can be used in an autonomous way and does not depend on CHAINS. It does not need any files outside the ones present in its own directory. It can be extracted and used to launch independent QOCT-RA calculations, and can be adapted to deal with other similar problematics. As such, this specific part of the documentation only explains how ``CONTROL LAUNCHER`` works in and of itself, and its integration into CHAINS is explained :doc:`elsewhere <chains.control_profile>`.

.. _control_directory:

Directory structure
===================

The ``CONTROL LAUNCHER`` `directory <https://github.com/niacobel/CHAINS/tree/master/control_launcher>`_ has and must keep the following structure:

.. code-block:: text

    control_launcher/
      ├── control_launcher.py
      ├── source_parser.py
      ├── transition_fcts.py
      ├── control_renderer.py
      ├── control_errors.py
      ├── clusters.yml
      └── templates/

As for what each file does, everything will be explained in more details in the other sections of this documentation, but here is a short summary:

- ``control_launcher.py`` is the main script itself, the one that needs to be executed.
- ``source_parser.py`` is the library of functions that define how to parse the source file, i.e. how to fetch from it the needed molecular properties.
- ``transition_fcts.py`` is the library of functions that define which transitions will be covered, by creating the various initial and target states files.
- ``control_renderer.py`` is the library of functions that define how to render the Jinja templates, i.e. how to create the parameters files and the job scripts.
- ``control_errors.py`` contains all the classes and functions defining :ref:`how to handle errors <control_errors>`.
- ``clusters.yml`` is the YAML file containing all the information specific to the different clusters, called the :ref:`clusters configuration file <control_clusters_file>`.
- ``templates`` is the directory containing all the Jinja templates that will be used by ``CONTROL LAUNCHER``. 

In CHAINS' case, the ``templates`` directory structure is:

.. code-block:: text

   └── templates/
         ├── param.nml.jinja
         └── qoctra_job.sh.jinja

where

- ``param.nml.jinja`` is the Jinja template for the QOCT-RA parameters file.
- ``qoctra_job.sh.jinja`` is the Jinja template for the QOCT-RA job script. It contains all the commands that will be executed through the job scheduler to perform the calculation.

.. _control_arguments:

Command line arguments
======================

.. argparse::
   :module: control_launcher
   :func: parser
   :prog: control_launcher.py
   :nodescription:

How does it work?
=================

The executable part of ``CONTROL LAUNCHER`` is the main script, ``control_launcher.py``. This is the one that must be called in the command line. The overall procedure follows three main steps: **system modelling**, **determining the transitions** and **jobs launching**. Those steps will be more thoroughly explained in the other sections of this documentation. As such, this subsection will only focus on the global procedure.

.. Important::

   An important file that will be often referenced throughout this documentation is the **YAML clusters configuration file** (``clusters.yml``). Rather than presenting it in its entirety at the beginning, the relevant bits of information will be introduced in the different sections, but you can have a full overview of that file in its :ref:`specific documentation <control_clusters_file>`.

Input files
-----------

There are two main input files for ``CONTROL LAUNCHER``:

- :guilabel:`-s / \\--source`, the **source file**.

   The file containing all the values for the different molecular properties needed by the control procedure, such as the energy of the states, the coupling elements between them, the transition dipole moments and so on.

- :guilabel:`-cf / \\--config`, the **configuration files**.

   The YAML files containing the parameters specific to your QOCT-RA calculation (number of time steps, number of iterations, etc.). Those files must have the .yml or .yaml extension. You can either indicate a specific file in the command line, or point towards a directory where there are multiple of those files. If you specify multiple configuration files, ``CONTROL LAUNCHER`` will process each transition-configuration combination. For example, if you have 4 possible transitions and 3 configuration files, you will end up with 12 launched jobs on your cluster.

Other arguments
---------------

There are three other required arguments for executing ``CONTROL LAUNCHER``:

- :guilabel:`-cl / \\--cluster_name`, the **name of the cluster** you are running on.

   This value must be the same as the one given in the :ref:`clusters configuration file <control_clusters_file>`, so that ``CONTROL LAUNCHER`` knows what you are referring to. *(This is case-sensitive!)*

.. Tip::

   This argument does not need to be the same name as the actual name of your machine. It is just a label used by ``CONTROL LAUNCHER`` to know which information to get from its clusters configuration file.

- :guilabel:`-p / \\--profile`, the **name of the profile** you want to run jobs with.

   The profile is a label used by ``CONTROL LAUNCHER`` to know which information to get from its different files. It specifies what type of parsing and rendering will be performed as well as the transitions that will be covered. As will be shown throughout this documentation, the profiles are defined in the :ref:`clusters configuration file <control_clusters_file>` and the name given in the command line must match one of the names given in the file. *(This is case-sensitive!)*.

- :guilabel:`-o / \\--out_dir`, the **output directory** 

   This is the jobs root directory, where each molecule directory will be created. Those directories will contain the data files and the job subdirectories, from where the jobs will be submitted to the job scheduler. See :ref:`control_out_dir_struct` for details.

There are also some optional arguments that can be used to adapt to some specific situations. Their description in the :ref:`command line arguments <control_arguments>` subsection should be self-explanatory.

.. _system_modelling:

First step: System modelling
----------------------------

``CONTROL LAUNCHER`` begins by generating some data files that will be used by QOCT-RA to model the system that needs to be controlled. This step involves parsing the source file to extract all the needed values from it, then manipulating those values to build an effective Hamiltonian describing the molecule. 

For more details on how this modelling process is performed, consult the :doc:`control_launcher.modelling` specific documentation.

.. _determining_transitions:

Second step: Determining the transitions
----------------------------------------

Once the system has been modelled, it is time to determinate the transitions that will be covered by the control procedure. For this, some data files still need to be created: the transition files, i.e. the initial states and the target states. 

For more details on how this process is performed, consult the :doc:`control_launcher.transitions` specific documentation.

Third step: Jobs launching
--------------------------

After having created all the data files needed to model the system, ``CONTROL LAUNCHER`` creates the job directories and files, then launches the actual jobs themselves. This step can be divided into three smaller steps:

- **Job scaling**: Using the number of states as a reference, the job scale for that molecule is evaluated, which will specify the calculation requirements accordingly (walltime, memory, etc.). For more details, consult the :doc:`control_launcher.job_scale` specific documentation.
- **Rendering the templates**: based on user-defined Jinja templates, ``CONTROL LAUNCHER`` creates the input files and the job script associated with our calculation. The content of those files is based on the information from the configuration files. For more details, consult the :doc:`control_launcher.rendering` specific documentation.
- **Submitting the job**: Now that everything has been prepared for the job, ``CONTROL LAUNCHER`` submits it to the job scheduler.

The exact command that will be executed for submitting the job is:

.. code-block:: console

    $ <submit_command> <delay_command> <job script>

where

- ``<submit_command>`` is the command which submits jobs to your job scheduler. In SLURM's case, it is the ``sbatch`` command. This must be indicated in the :ref:`clusters configuration file <control_clusters_file>`: 

   .. code-block:: yaml

      mycluster:
        submit_command: <submit_command>

   where ``mycluster`` is the name of your cluster (the same that was given as the :guilabel:`-cl / \\--cluster_name` command line argument).

- ``<delay_command>`` is an optional command that can delay the submission of a particular job, which can prove useful if you want to prioritize certain job sizes (consult the :doc:`control_launcher.job_scale` specific documentation for details). In SLURM's case, this is covered by the ``--begin`` argument.
- ``<job script>`` is the name of the file that will be created through the :doc:`rendering process <control_launcher.rendering>`. It contains the commands needed by the job scheduler to run the calculation on the cluster.

For example, if we want to run a QOCT-RA calculation on a SLURM cluster, but delay the submission of this job by 60 seconds, the command executed by ``CONTROL LAUNCHER`` might look like:

.. code-block:: console

    $ sbatch --begin=now+60 qoctra_job.sh

Once the job has been submitted, ``CONTROL LAUNCHER`` will proceed to the next configuration file with the same transition. Once all the configuration files have been treated, it will proceed to the next transition and treat again all the configuration files for that transition. At the end of the execution, barring any problems, a job will have been launched for each transition-configuration combination.

.. _control_out_dir_struct:

Output directory structure
--------------------------

If we have for example 2 transitions and 2 configuration files, once the execution of ``CONTROL LAUNCHER`` has ended, the structure of the output directory (given as the :guilabel:`-o / \\--out_dir` command line argument) might look like:

.. code-block:: text

    out_dir/ 
      └── source/
            └── data/
                  ├── source_file
                  ├── source.log                
                  └── all the data files
            └── transition1_config1/
                  ├── config1.yml
                  ├── transition1_config1.log
                  ├── job_script.sh
                  └── param.nml
            └── transition1_config2/
                  ├── config2.yml
                  ├── transition1_config2.log
                  ├── job_script.sh
                  └── param.nml
            └── transition2_config1/
                  ├── config1.yml
                  ├── transition2_config1.log
                  ├── job_script.sh
                  └── param.nml
            └── transition2_config2/
                  ├── config2.yml
                  ├── transition2_config2.log
                  ├── job_script.sh
                  └── param.nml

where 

- ``source`` is the directory created by ``CONTROL LAUNCHER`` and named after the source file (minus a possible extension)
- ``data`` is the directory containing all the data files created during the :ref:`system modelling <system_modelling>` and :ref:`determining the transitions <determining_transitions>` steps.
- ``source_file`` is a copy of the source file.
- ``source.log`` is an output file containing the details of the treatment of this source file by ``CONTROL LAUNCHER`` (the extracted molecular properties, the considered transitions, etc.).
- ``transitionX_configX`` is the job subdirectory from which the job will be submitted to the job scheduler.
- ``configX.yml`` is a copy of the configuration file.
- ``job_script.sh`` and ``param.nml`` are the files created by the :doc:`rendering process <control_launcher.rendering>`.
- ``transitionX_configX.log`` is an output file containing the details of the treatment of this transition-configuration combination by ``CONTROL LAUNCHER`` (the used job scale, the files created, etc.).

.. Hyperlink targets

.. _QOCT-RA: https://gitlab.com/dynaq.cqp/QOCT-RA