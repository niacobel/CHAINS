*********************
Other important files
*********************

.. _clusters_file:

Clusters configuration file
===========================

The clusters configuration file is a YAML file containing all the information about the clusters, such as the command used to submit jobs, the location of the programs you want to run, the command to run them, the scaling function and the job scales. If you don't know what those parameters are, they have been introduced in the previous sections of this documentation (see :ref:`submitting step <submitting_step>`, :doc:`job scaling <abin_launcher.job_scale>` and :doc:`rendering the templates <abin_launcher.rendering>`). Here is a full overview of what the file might look like:

.. code-block:: yaml

   myclusterA:
     subcommand: value
     progs:
       myprog1:
         set_env: value
         command: value
         scaling_function: name-of-scaling-function
         job_scales:
           - label: scale1
             scale_limit: value
             partition_name: value
             time: value
             cores: value
             mem_per_cpu: value
             delay_command: value
           - label: scale2
             scale_limit: value
             partition_name: value
             time: value
             cores: value
             mem_per_cpu: value
             delay_command: value
           - ...
       myprog2:
         set_env: value
         command: value
         scaling_function: name-of-scaling-function
         job_scales:
           - label: scale1
             scale_limit: value
             partition_name: value
             time: value
             cores: value
             mem_per_cpu: value
             delay_command: value
           - label: scale2
             scale_limit: value
             partition_name: value
             time: value
             cores: value
             mem_per_cpu: value
             delay_command: value
           - ...

   myclusterB:
     ...

If you want a more concrete example, let's consider the following situation:

- Two clusters who use SLURM as the job scheduler, named slurm_cluster_A and slurm_cluster_B
- We want to run two programs: ORCA and Q-CHEM
- Both programs use ``total_nb_elec`` as the scaling function
- ORCA is available on both clusters, but Q-CHEM is specific to the slurm_cluster_B

This is what the file might look like in this situation:

.. code-block:: yaml

   slurm_cluster_A:
     subcommand: sbatch
     progs:
       orca:
         set_env: module load ORCA/4.1.0-OpenMPI-3.1.3
         command: /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca
         scaling_function: total_nb_elec
         job_scales:
           - label: tiny
             scale_limit: 50
             partition_name: batch
             time: 0-00:10:00
             cores: 4 
             mem_per_cpu: 500 # in MB
             delay_command:
           - label: small
             scale_limit: 500
             partition_name: batch
             time: 1-00:00:00
             cores: 8
             mem_per_cpu: 500 # in MB
             delay_command:
           - label: medium
             scale_limit: 1000
             partition_name: batch
             time: 2-00:00:00
             cores: 8
             mem_per_cpu: 2000 # in MB
             delay_command: --begin=now+60
           - label: big
             scale_limit: 1500
             partition_name: batch
             time: 3-00:00:00
             cores: 16
             mem_per_cpu: 4000 # in MB
             delay_command: --begin=now+120

   slurm_cluster_B:
     subcommand: sbatch
     progs:
       orca:
         set_env: module load orca/4.0.1.2
         command: /usr/local/orca/orca_4_0_1_2_linux_x86-64_openmpi202/orca
         scaling_function: total_nb_elec
         job_scales:
           - label: tiny
             scale_limit: 50
             partition_name: default
             time: 0-00:20:00
             cores: 4 
             mem_per_cpu: 500 # in MB
             delay_command: 
           - label: small
             scale_limit: 500
             partition_name: default
             time: 1-10:00:00
             cores: 8
             mem_per_cpu: 500 # in MB
             delay_command: 
           - label: medium
             scale_limit: 1000
             partition_name: default
             time: 3-00:00:00
             cores: 8
             mem_per_cpu: 2000 # in MB
             delay_command: --begin=now+60
       qchem:
         set_env: module load Q-Chem/5.3.0-SHMEM
         command: srun qchem
         scaling_function: total_nb_elec
         job_scales:
           - label: tiny
             scale_limit: 100
             partition_name: default
             time: 0-00:20:00
             cores: 4 
             mem_per_cpu: 500 # in MB
             delay_command:
           - label: small
             scale_limit: 750
             partition_name: default
             time: 1-00:00:00
             cores: 8
             mem_per_cpu: 1000 # in MB
             delay_command:
           - label: medium
             scale_limit: 1500
             partition_name: default
             time: 3-00:00:00
             cores: 8
             mem_per_cpu: 2000 # in MB
             delay_command: --begin=now+60
           - label: big
             scale_limit: 2000
             partition_name: long
             time: 8-00:00:00
             cores: 16
             mem_per_cpu: 4000 # in MB
             delay_command: --begin=now+120

This is the bare minimum required by ``ABIN LAUNCHER``, but you can add as many keys as you want depending on your needs. Just remember to adjust the ``scale_limit`` of your job scales if you change your scaling function. Otherwise, those numbers won't make sense.

.. _abin_errors:

Errors handling
===============

When adding a :ref:`rendering function <rendering_fct>` or another custom function to ``ABIN LAUNCHER``, having a way to handle errors is definitely useful. In ``ABIN LAUNCHER``, this is managed by the ``abin_errors.py`` file. It is somewhat basic but should be enough to cover your needs.

Custom exception
----------------

A custom exception class has been created to handle errors specific to ``ABIN LAUNCHER``, in the ``abin_errors.py`` file:

.. autoclass:: abin_errors.AbinError 

Feel free to call it when you want to prevent predictable errors from happening (missing file, incorrect value, etc.) by simply using

.. code-block:: python

   raise abin_errors.AbinError (<message>)

where ``<message>`` is a proper error message describing the error. Those raised exceptions wil be handled ``ABIN LAUNCHER``, which will then either abort the execution or skip the incriminated geometry or configuration file, depending on where the error occurred.

Checking the existence of files and directories
-----------------------------------------------

In order to easily check if specific files or directories exist, a function has been defined in the ``abin_errors.py`` file:

.. autofunction:: abin_errors.check_abspath

Let's say we want to check if our periodic table file is still there, we can use the code:

.. code-block:: python

  mendeleev_file = abin_errors.check_abspath(os.path.join(code_dir,"mendeleev.yml"),"Mendeleev periodic table YAML file","file")

This line of code will check if there is a file named ``mendeleev.yml`` in ``ABIN LAUNCHER``'s directory (``code_dir``) and if it is indeed a file (and not a directory for example). If there is, it will return the absolute path towards that file (useful for referencing that file later in the script no matter where the current directory is). Otherwise, it will raise an exception and specify the context as "Mendeleev periodic table YAML file" for easy tracking, which might result in the following error message:

.. code-block:: text

   Something went wrong when checking the path  ~/CHAINS/abin_launcher/mendeleev.yml
   Context:  Mendeleev periodic table YAML file
   ERROR: ~/CHAINS/abin_launcher/mendeleev.yml does not seem to exist.
