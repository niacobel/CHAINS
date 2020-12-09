***********
Job scaling
***********

What is job scaling and why does it matter?
===========================================

"Job scaling" consists here to automatically match the computing resources requirement (time limit, number of CPUs, memory, etc.) to the size of the job. Without job scaling, you might either waste resources or not allocate enough of them. In the former scenario, your user fairshare will become uselessly big, and your subsequent jobs might wait a long time before starting. In the latter scenario, your job might just simply crash.

How does it work?
=================

This process first assigns a value to the molecule that will reflect the job size and complexity. That value is what we call the **scale index**. The way ``ABIN LAUNCHER`` computes this value is by calling a function defined in the ``scaling_fcts.py`` file, called a **scaling function**. Then, this scale index will be compared to what we defined as **job scales**. Those scales are a set of computing resources parameters associated with a range of scale index values. 

For example, let's say we want to run a geometry optimization on the |Si36Ge11H60| molecule. The scale index has been computed as 916 and the job scales are defined as follows:

.. code-block:: yaml

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
        time: 1-00:00:00
        cores: 8
        mem_per_cpu: 500 # in MB
      - 
        label: medium
        scale_limit: 1000
        time: 2-00:00:00
        cores: 8
        mem_per_cpu: 2000 # in MB
      - 
        label: big
        scale_limit: 1500
        time: 3-00:00:00
        cores: 16
        mem_per_cpu: 4000 # in MB

The ``scale_limit`` key defines the upper limit of that job scale for the scale index. This means our molecule is too big for the ``tiny`` and ``small`` scales, which have an upper limit of 50 and 500, respectively. It will then uses the resources defined in the ``medium`` scale, which are: a time limit of 2 days, 8 CPUs and 2000 MB of memory per CPU.

Obviously, the key part of this process lies in the quality of the job scales definition. The finer they are, the better the scaling will be. Since this is highly dependent on the program you want to run and the cluster on which it will be running, you will need to do extensive testing on your part. If your cluster uses SLURM as a job scheduler, you might want to take a look at the :doc:`abin_launcher.benchmark` section.

.. |Si36Ge11H60| replace:: Si\ :sub:`36`\ Ge\ :sub:`11`\ H\ :sub:`60`\ 

.. _scaling_fcts:

Scaling functions
=================

.. _scaling_fcts_definition:

General definition
------------------

All the scaling functions are defined in the ``scaling_fcts.py`` file and obey some restrictions, in order to be callable by ``ABIN LAUNCHER``:

- They only take two dictionaries as arguments: the content of the ``mendeleev.yml`` file and the ``file_data`` variable, as built by the :doc:`scanning function <abin_launcher.scan>`.
- They return an integer or a float, that will act as the scale index.
  
If a problem arises when computing the scale index, an ``AbinError`` exception is raised with a proper error message (see :ref:`how to handle errors <abin_errors>` for more details).

Choosing a scaling function
---------------------------

The scaling function that will be called by ``ABIN LAUNCHER`` is the one associated with the ``scaling_function`` YAML key defined in the :ref:`clusters configuration file <abin_clusters_file>`:

.. code-block:: yaml

   mycluster:
      profiles:
         myprofile1:
            scaling_function: name-of-scaling-function
         myprofile2:
            scaling_function: name-of-scaling-function

where ``mycluster`` corresponds to the name of your cluster (given as a :ref:`command line argument <abin_arguments>`) while ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to run (such as ``orca`` or ``qchem``). This way, a different scaling function can be assigned to each profile.

Total number of electrons
-------------------------

.. autofunction:: scaling_fcts.total_nb_elec

This function starts by defining a sub-function called ``get_nb_elec_for_element`` that makes use of Mendeleev's Periodic Table (``mendeleev.yml``) in order to get the atomic number of an atom type. Then, the main function gets the different atom types and their respective amount from the ``chemical_formula`` key in ``file_data``. For each atom type, the function multiplies the associated number of electrons (obtained through the sub-function) with the amount of electrons of that type in the molecule. It finally sums up the value obtained for each atom type and returns it as the scale index.

Example for the |C3H8| molecule:

.. code-block:: text

   ---------------------------------------------------------------------
   Atom Type    Atomic Number    Number of atoms     Number of electrons   
   ---------------------------------------------------------------------
   H            1                8                   8                     
   C            6                3                   18                    
   ---------------------------------------------------------------------
   Total                         11                  26                    
   ---------------------------------------------------------------------
   
   Scale index:                                      26    

.. |C3H8| replace:: C\ :sub:`3`\ H\ :sub:`8`\ 

Total number of atoms
---------------------

.. autofunction:: scaling_fcts.total_nb_atoms

This function simply sums up all the values of the keys in the ``chemical_formula`` key in ``file_data``. It does not make use of ``mendeleev.yml`` but it is still passed as an argument to satisfy the calling restrictions of the scaling functions (see :ref:`scaling_fcts_definition`).

.. _job_scales:

Job scales
==========

The job scales must be defined as follows in the ``job_scales`` key in the :ref:`clusters configuration file <abin_clusters_file>`:

.. code-block:: yaml

   myclusterA:
     profiles:
       myprofile1:
         scaling_function: name-of-scaling-function
         job_scales:     
           - 
             label: scale1
             scale_limit: value
             time: value
             cores: value
             mem_per_cpu: value
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
         scaling_function: name-of-scaling-function
         job_scales:     
           - 
             label: scale1
             scale_limit: value
             time: value
             cores: value
             mem_per_cpu: value
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
     profiles:
       myprofile1:
         scaling_function: name-of-scaling-function
         job_scales:     
           - 
             label: scale1
             scale_limit: value
             time: value
             cores: value
             mem_per_cpu: value
             partition_name: value  # This is optional
           - 
             label: scale2
             scale_limit: value
             time: value
             cores: value
             mem_per_cpu: value
             delay_command: value   # This is optional
           - 
             ...
  
       myprofile2:
         scaling_function: name-of-scaling-function
         job_scales:     
           - 
             label: scale2
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
           - 
             ...

where

- ``myclusterA`` and ``myclusterB`` are the names of your clusters (given as a :ref:`command line argument <abin_arguments>`). This way, different job scales can be assigned to each cluster.
- ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to run (such as ``orca`` or ``qchem``, given as a :ref:`command line argument <abin_arguments>`). This way, different job scales can be assigned to each profile.
- ``label``, ``scale_limit``, ``time``, ``cores`` and ``mem_per_cpu`` are all **mandatory keys**, specifying the resources requirements of the jobs.
- ``partition_name`` is an optional key containing the name of the cluster partition on which the job will be running.
- ``delay_command`` is an optional key that lets you delay the submission of the jobs. For example, by delaying the bigger jobs, you can prioritize the launch of small calculations first. On SLURM, this is handled by the ``--begin`` argument of the ``sbatch`` command, see here_.

You can have as many job scales as you want, and they don't need to be defined in ascending order of scale index limits. ``ABIN LAUNCHER`` will automatically sort them before starting to scan the geometry files. Just remember to adjust the ``scale_limit`` of your job scales if you change your scaling function. Otherwise, those numbers won't make sense.

.. Hyperlink targets

.. _here: https://slurm.schedmd.com/sbatch.html