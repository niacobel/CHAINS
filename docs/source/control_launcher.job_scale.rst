***********
Job scaling
***********

What is job scaling and why does it matter?
===========================================

"Job scaling" consists here to automatically match the computing resources requirement (time limit, memory, etc.) to the size of the job. Without job scaling, you might either waste resources or not allocate enough of them. In the former scenario, your user fairshare will become uselessly big, and your subsequent jobs might wait a long time before starting. In the latter scenario, your job might just simply crash.

How does it work?
=================

This process first assigns a value to the molecule that will reflect the job size and complexity. That value is what we call the **scale index**. For ``CONTROL LAUNCHER``, this value is the number of states involved in the control procedure. Then, this scale index will be compared to what we defined as **job scales**. Those scales are a set of computing resources parameters associated with a range of scale index values. 

For example, let's consider that we are working with 24 electronic states of the |Si17H36| molecule: the ground state, twelve excited triplet states and twelve excited singlet states. The scale index is then 24 and the job scales are defined as follows:

.. code-block:: yaml

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

The ``scale_limit`` key defines the upper limit of that job scale for the scale index. This means our number of states is too big for the ``small`` scale, which has an upper limit of 20. It will then uses the resources defined in the ``medium`` scale, which are: a time limit of 2 days and 2500 MB of memory.

Obviously, the key part of this process lies in the quality of the job scales definition. The finer they are, the better the scaling will be. Since this is highly dependent on the program you want to run and the cluster on which it will be running, you will need to do extensive testing on your part.

.. |Si17H36| replace:: Si\ :sub:`17`\ H\ :sub:`36`\ 

.. _control_job_scales:

Job scales
==========

The job scales must be defined as follows in the ``job_scales`` key in the :ref:`clusters configuration file <control_clusters_file>`:

.. code-block:: yaml

   myclusterA:
     profiles:
       myprofile1:
         job_scales:     
           - 
             label: scale1
             scale_limit: value
             time: value
             memory: value
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
         job_scales:     
           - 
             label: scale1
             scale_limit: value
             time: value
             memory: value
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
     profiles:
       myprofile1:
         job_scales:     
           - 
             label: scale1
             scale_limit: value
             time: value
             memory: value
             partition_name: value  # This is optional
           - 
             label: scale2
             scale_limit: value
             time: value
             memory: value
             delay_command: value   # This is optional
           - 
             ...
  
       myprofile2:
         job_scales:     
           - 
             label: scale2
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
           - 
             ...

where

- ``myclusterA`` and ``myclusterB`` are the names of your clusters (given as a :ref:`command line argument <control_arguments>`). This way, different job scales can be assigned to each cluster.
- ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to use (such as ``chains_qoctra`` or ``basic_qoctra``, given as a :ref:`command line argument <control_arguments>`). This way, different job scales can be assigned to each profile.
- ``label``, ``scale_limit``, ``time`` and ``memory`` are all **mandatory keys**, specifying the resources requirements of the jobs.
- ``partition_name`` is an optional key containing the name of the cluster partition on which the job will be running.
- ``delay_command`` is an optional key that lets you delay the submission of the jobs. For example, by delaying the bigger jobs, you can prioritize the launch of small calculations first. On SLURM, this is handled by the ``--begin`` argument of the ``sbatch`` command, see here_.

You can have as many job scales as you want, and they don't need to be defined in ascending order of scale index limits. ``CONTROL LAUNCHER`` will automatically sort them before starting to parse the source file.

.. Hyperlink targets

.. _here: https://slurm.schedmd.com/sbatch.html