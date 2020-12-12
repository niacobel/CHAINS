*********************
Other important files
*********************

.. _control_clusters_file:

Clusters configuration file
===========================

The clusters configuration file, named ``clusters.yml``, is a YAML file where all the information about the clusters is specified. This includes the command used to submit jobs on the cluster as well as all the profiles. Those profiles define the parsing function, the transition function, the rendering function, the location of the compilers needed to run QOCT-RA and the job scales. If you don't know what those parameters are, they have been introduced in the previous sections of this documentation (see :ref:`jobs launching step <jobs_launching>`, :doc:`system modelling <control_launcher.modelling>`, :doc:`determining the transitions <control_launcher.transitions>`, :doc:`job scaling <control_launcher.job_scale>` and :doc:`rendering the templates <control_launcher.rendering>`). Here is what the structure of this file looks like:

.. code-block:: yaml

   myclusterA:
     submit_command: value
     profiles:
       myprofile1:
         parsing_function: name-of-parsing-function
         transition_function: name-of-transition-function
         rendering_function: name-of-rendering-function
         set_env: value
         job_scales:
           - 
             label: scale1
             scale_limit: value
             time: value
             memory: value
             partition_name: value  # This is optional
             delay_command: value   # This is optional
           - 
             label: scale2
             scale_limit: value
             time: value
             memory: value
             partition_name: value  # This is optional
             delay_command: value   # This is optional
           - 
             ...            
       myprofile2:
         parsing_function: name-of-parsing-function
         transition_function: name-of-transition-function
         rendering_function: name-of-rendering-function
         set_env: value
         job_scales:
           - 
             label: scale1
             scale_limit: value
             time: value
             memory: value
             partition_name: value  # This is optional
             delay_command: value   # This is optional
           - 
             label: scale2
             scale_limit: value
             time: value
             memory: value
             partition_name: value  # This is optional
             delay_command: value   # This is optional
           - 
             ...

   myclusterB:
     ...

where

- ``myclusterA`` and ``myclusterB`` are the names of your clusters (given as a :ref:`command line argument <abin_arguments>`).
- ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to run (also given as a :ref:`command line argument <abin_arguments>`).

If you want a more concrete example, let's consider the following situation:

- A cluster who use SLURM as the job scheduler, named ``lemaitre3``
- The ``qchem_gt_opm`` profile with ``qchem_tddft``, ``proj_ground_to_triplet`` and ``opm_render`` as the parsing, transition and rendering functions, respectively.

This is what the file might look like in this situation:

.. code-block:: yaml

   lemaitre3:
     submit_command: sbatch
     profiles:
       qchem_gt_opm:
         parsing_function: qchem_tddft
         transition_function: proj_ground_to_triplet
         rendering_function: opm_render
         set_env: module load OpenBLAS/0.3.7-GCC-8.3.0                  
         job_scales:
         - 
           label: small
           scale_limit: 20
           time: 1-00:00:00
           memory: 2000 # in MB
         - 
           label: medium
           scale_limit: 50
           time: 2-00:00:00
           memory: 2500 # in MB
         - 
           label: big
           scale_limit: 100
           time: 5-00:00:00
           memory: 3000 # in MB

This is what a basic example looks like, but you can add as many keys as you want, depending on your needs.

.. _control_errors:

Errors handling
===============

When adding a :ref:`rendering function <control_rendering_fct>` or another custom function to ``CONTROL LAUNCHER``, having a way to handle errors is definitely useful. In ``CONTROL LAUNCHER``, this is managed by the ``control_errors.py`` file. It is somewhat basic but should be enough to cover your needs.

Custom exception
----------------

A custom exception class has been created to handle errors specific to ``CONTROL LAUNCHER``, in the ``control_errors.py`` file:

.. autoclass:: control_errors.ControlError 

Feel free to raise it when you want to prevent predictable errors from happening (missing file, incorrect value, etc.) by simply using

.. code-block:: python

   raise control_errors.ControlError ("my message here")

Those raised exceptions wil be caught by ``CONTROL LAUNCHER``, which will then either abort the execution or skip the incriminated configuration file, depending on where the error occurred.

Checking the existence of files and directories
-----------------------------------------------

In order to easily check if specific files or directories exist, the ``check_abspath`` function has been defined in the ``control_errors.py`` file:

.. autofunction:: control_errors.check_abspath

As an example, let's say we want to check if our clusters configuration file is still there, we can use the code:

.. code-block:: python

  clusters_file = control_errors.check_abspath(os.path.join(code_dir,"clusters.yml"),"YAML clusters configuration file","file")

This will check if there is a file named ``clusters.yml`` in ``CONTROL LAUNCHER``'s directory (``code_dir``) and if it is indeed a file (and not a directory for example). 

- If there is, it will return the absolute path towards that file (useful for referencing that file later in the script, no matter where the current directory is). 
- Otherwise, it will raise an exception and specify the context as "YAML clusters configuration file" for easy tracking, which will result in an error message of the form:

.. code-block:: console

   Something went wrong when checking the path  ~/CHAINS/control_launcher/clusters.yml
   Context:  YAML clusters configuration file
   ERROR: ~/CHAINS/control_launcher/clusters.yml does not seem to exist.

.. Hyperlink targets
