*********************
Other important files
*********************

.. _abin_clusters_file:

Clusters configuration file
===========================

The clusters configuration file, named ``clusters.yml``, is a YAML file where all the information about the clusters is specified. This includes the command used to submit jobs on the cluster as well as all the profiles. Those profiles define the rendering function, the location of the programs you want to run, the command to run them, the scaling function and the job scales. If you don't know what those parameters are, they have been introduced in the previous sections of this documentation (see :ref:`submitting step <submitting_step>`, :doc:`job scaling <abin_launcher.job_scale>` and :doc:`rendering the templates <abin_launcher.rendering>`). Here is what the structure of this file looks like:

.. code-block:: yaml

   myclusterA:
     submit_command: value
     profiles:
       myprofile1:
         rendering_function: name-of-rendering-function
         set_env: value
         command: value
         scaling_function: name-of-scaling-function
         job_scales:
           - 
             label: scale1
             scale_limit: value
             time: value
             cores: value
             mem_per_cpu: value
             partition_name: value  # This is optional
             delay_command: value   # This is optional
           - 
             label: scale2
             scale_limit: value
             time: value
             cores: value
             mem_per_cpu: value
             partition_name: value  # This is optional
             delay_command: value   # This is optional
           - 
             ...
       myprofile2:
         rendering_function: name-of-rendering-function
         set_env: value
         command: value
         scaling_function: name-of-scaling-function
         job_scales:
           - 
             label: scale1
             scale_limit: value
             time: value
             cores: value
             mem_per_cpu: value
             partition_name: value  # This is optional
             delay_command: value   # This is optional
           - 
             label: scale2
             scale_limit: value
             time: value
             cores: value
             mem_per_cpu: value
             partition_name: value  # This is optional
             delay_command: value   # This is optional
           - 
             ...

   myclusterB:
     ...

where

- ``myclusterA`` and ``myclusterB`` are the names of your clusters (given as a :ref:`command line argument <abin_arguments>`).
- ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to use (also given as a :ref:`command line argument <abin_arguments>`).

If you want a more concrete example, let's consider the following situation:

- Two clusters who use SLURM as the job scheduler, named ``vega`` and ``lemaitre3``
- Two profiles for ``vega``: ``orca`` and ``qchem`` with ``orca_render`` and ``qchem_render`` as rendering functions, respectively.
- One profile for ``lemaitre3`` who also uses the ``orca`` profile but with different commands to load and execute ORCA_, as well as different job scales.
- All the profiles use ``total_nb_elec`` as the scaling function.

This is what the file might look like in this situation:

.. code-block:: yaml

   vega:
     submit_command: sbatch
     profiles:
       orca:
         rendering_function: orca_render
         set_env: module load ORCA/4.0.0.2-OpenMPI-2.0.2
         command: /apps/brussel/interlagos/software/ORCA/4.0.0.2-OpenMPI-2.0.2/orca
         scaling_function: total_nb_elec
         job_scales:
           - 
             label: tiny
             scale_limit: 50
             time: 0-00:20:00
             cores: 4 
             mem_per_cpu: 500 # in MB 
           - 
             label: small
             scale_limit: 500
             time: 1-10:00:00
             cores: 8
             mem_per_cpu: 500 # in MB 
           - 
             label: medium
             scale_limit: 1000
             time: 3-00:00:00
             cores: 8
             mem_per_cpu: 2000 # in MB
             delay_command: --begin=now+60
       qchem:
         rendering_function: qchem_render
         set_env: module load Q-Chem-5.2.1-intel-2019b-mpich3
         command: srun qchem
         scaling_function: total_nb_elec
         job_scales:
           - 
             label: tiny
             scale_limit: 100
             time: 0-00:20:00
             cores: 4 
             mem_per_cpu: 500 # in MB
           - 
             label: small
             scale_limit: 750
             time: 1-00:00:00
             cores: 8
             mem_per_cpu: 1000 # in MB
           - 
             label: medium
             scale_limit: 1500
             time: 3-00:00:00
             cores: 8
             mem_per_cpu: 2000 # in MB
             delay_command: --begin=now+60
           - 
             label: big
             scale_limit: 2000
             partition_name: long
             time: 8-00:00:00
             cores: 16
             mem_per_cpu: 4000 # in MB
             delay_command: --begin=now+120

   lemaitre3:
     submit_command: sbatch
     profiles:
       orca:
         rendering_function: orca_render
         set_env: module load ORCA/4.1.0-OpenMPI-3.1.3
         command: /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca
         scaling_function: total_nb_elec
         job_scales:
           - 
             label: tiny
             scale_limit: 50
             time: 0-00:10:00
             cores: 4 
             mem_per_cpu: 500 # in MB
           - 
             label: small
             scale_limit: 500
             partition_name: batch
             time: 1-00:00:00
             cores: 8
             mem_per_cpu: 500 # in MB
           - 
             label: medium
             scale_limit: 1000
             partition_name: batch
             time: 2-00:00:00
             cores: 8
             mem_per_cpu: 2000 # in MB
             delay_command: --begin=now+60
           - 
             label: big
             scale_limit: 1500
             partition_name: batch
             time: 3-00:00:00
             cores: 16
             mem_per_cpu: 4000 # in MB
             delay_command: --begin=now+120

This is what a basic example looks like, but you can add as many keys as you want, depending on your needs.

.. _abin_errors:

Errors handling
===============

When adding a :ref:`rendering function <rendering_fct>` or another custom function to ``ABIN LAUNCHER``, having a way to handle errors is definitely useful. In ``ABIN LAUNCHER``, this is managed by the ``abin_errors.py`` file. It is somewhat basic but should be enough to cover your needs.

Custom exception
----------------

A custom exception class has been created to handle errors specific to ``ABIN LAUNCHER``, in the ``abin_errors.py`` file:

.. autoclass:: abin_errors.AbinError 

Feel free to raise it when you want to prevent predictable errors from happening (missing file, incorrect value, etc.) by simply using

.. code-block:: python

   raise abin_errors.AbinError ("my message here")

Those raised exceptions wil be caught by ``ABIN LAUNCHER``, which will then either abort the execution or skip the incriminated geometry or configuration file, depending on where the error occurred.

Checking the existence of files and directories
-----------------------------------------------

In order to easily check if specific files or directories exist, the ``check_abspath`` function has been defined in the ``abin_errors.py`` file:

.. autofunction:: abin_errors.check_abspath

As an example, let's say we want to check if our periodic table file is still there, we can use the code:

.. code-block:: python

  mendeleev_file = abin_errors.check_abspath(os.path.join(code_dir,"mendeleev.yml"),"Mendeleev periodic table YAML file","file")

This will check if there is a file named ``mendeleev.yml`` in ``ABIN LAUNCHER``'s directory (``code_dir``) and if it is indeed a file (and not a directory for example). 

- If there is, it will return the absolute path towards that file (useful for referencing that file later in the script, no matter where the current directory is). 
- Otherwise, it will raise an exception and specify the context as "Mendeleev periodic table YAML file" for easy tracking, which will result in an error message of the form:

.. code-block:: console

   Something went wrong when checking the path  ~/CHAINS/abin_launcher/mendeleev.yml
   Context:  Mendeleev periodic table YAML file
   ERROR: ~/CHAINS/abin_launcher/mendeleev.yml does not seem to exist.

.. Hyperlink targets

.. _ORCA: https://www.faccts.de/orca/
.. _Q-CHEM: https://www.q-chem.com/