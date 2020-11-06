*******************
Benchmarking option
*******************

.. warning::

   This section is still a work in progress and some information might be missing or temporarily erroneous.

Why is it useful?
=================

As mentioned in the :doc:`abin_launcher.job_scale` section, one of the key parts of the scaling process is the definition of your job scales. A good amount of finely tuned job scales will allow you to avoid wasting resources, thus diminishing your fairshare usage of the cluster, making your jobs spend less time waiting in the queue. Since this is highly dependent on the program you want to run and the cluster on which it will be running, you will need to do extensive testing on your part. 

With that said, if your cluster uses SLURM as a job scheduler, we might be able to offer some help in that regard through a benchmarking option. This option will allow us to see the efficiency of our resources requirements, by calculating how much of the allocated resources were actually used by the job.

Important SLURM commands
========================

The ``seff`` command: a possible alternative
--------------------------------------------

The first option is the ``seff`` command (see `source code <https://github.com/SchedMD/slurm/blob/master/contribs/seff/seff>`_ and `example <https://sites.google.com/a/case.edu/hpcc/jobs/slurm-command-overview/seff>`_). While it pretty much does exactly what we want, there are some problems with it:

- It requires SLURM 15.08 or a more recent version, which might cause problems with some older machines.
- It does not give reliable statistics when the job is running, so a crontab task or similar* will be needed in order to automatically check the resources usage after the job has ended.
- It computes the CPU and memory efficiencies, but nothing is done about the time efficiency.

Other than that, it still is a good option so feel free to use it if it satisfies your needs!

\* If your cluster administrator allows it, you can also use the ``strigger`` command to trigger a script executing the ``seff`` command at the end of the job, consult the `documentation <https://slurm.schedmd.com/strigger.html>`_ for details.

Our choice: the ``sacct`` and ``squeue`` commands
-------------------------------------------------

Given ``seff``'s limitations, we decided to go with two other commands here, ``sacct`` and ``squeue``. The ``sacct`` command can directly fetch job accounting data from the SLURM database. Depending on your cluster configuration however, reaching the SLURM database can necessitate an Internet access, which the computing nodes might not have. This command will then only be used from the login nodes, when the job has finished. On the other hand, the ``squeue`` command gets its information locally, which means we can use it for the information that will be obtained through the job instructions file, while the job is still running.

Every data about our jobs will be obtained through those two commands. If you want to customize the benchmarking procedure, you should familiarize yourself with them on SLURM's official documentation (`sacct docs <https://slurm.schedmd.com/sacct.html>`_ and `squeue docs <https://slurm.schedmd.com/squeue.html>`_). 

Required files
==============

This process requires three files:

- A Jinja template, named ``benchmark.jinja``, placed in the ``templates`` directory of ``ABIN LAUNCHER``. It is an extension of the job instructions template.
- A Python script, named ``benchmark.py``, which must be placed in ``ABIN LAUNCHER``'s directory.
- A Shell script, named ``cron_benchmark.sh``, which must also be placed in ``ABIN LAUNCHER``'s directory. It will be executed through a crontab task and make the link between the first two files.

How does it work?
=================

The role of the Jinja template
------------------------------

At the end of the job instructions, some additional commands, provided by the ``benchmark.jinja`` template, will write a new line in a **temporary CSV file**, containing the following information:

- The name of the program and the name of the cluster
- The chosen job scale and its associated resources requirements
- The job ID and name
- The scaling function and the computed scale index
- The number of nodes and the nodes list
- Four dates: the submit date, the eligible date, the start date and the end date (which is the current date at which the job ends since those instructions are executed at the end of the job script)

Some of those information are provided directly by ``ABIN LAUNCHER`` while the others are obtained through the ``squeue`` command.

The role of the crontab script
------------------------------

Meanwhile, on the cluster, the ``cron_benchmark.sh`` Shell script will be periodically executed through a crontab task to check the existence of that temporary CSV file. If the file exists, the script will archive it then execute the ``benchmark.py`` Python script to scan and process its content.

The role of the Python script
-----------------------------

The ``benchmark.py`` Python script will write a new line in a **final CSV file**. That line is a repeat of what was already present in the temporary file, with the following additions:

- The "reserved", or queued, time (time between the submit date and the eligible date)
- The elapsed time (duration of the job)
- The **time efficiency** (percentage of used time vs required time)
- The maximum used memory (in MB)
- The **memory efficiency** (percentage of used memory vs required memory)
- The total time used by all of the CPUs
- The "wall CPU", which is the total amount of time CPUs could have used (derived from the duration of the job and the number of CPUs)
- The **CPU efficiency** (percentage of total time used by the CPUs vs "wall CPU")

The only thing left to do is then to make a copy of that final CSV file on your local computer and open it in your favorite spreadsheet like Microsoft Excel!

How do I use it?
================

Prepare the Jinja template
--------------------------

First of all, make sure the ``benchmark.jinja`` template is present in the ``templates`` directory of ``ABIN LAUNCHER``. Then add the following line at the end of your job instructions template (which should be in the same directory):

.. code-block:: jinja

   {% include "benchmark.jinja" %}

Since that template requires some specific variables, add the following code to your :ref:`rendering function <rendering_fct>` *after* having defined your ``render_vars`` dictionary for your job instructions file, but *before* calling the ``jinja_render`` function for that file:

.. code-block:: python

    render_vars.update({
        "benchmark_dir" : "path/to/benchmark_dir",
        "prog" : job_specs['prog'],
        "cluster_name" : job_specs['cluster_name'],
        "jobscale_label" : job_specs['scale_label'],
        "job_walltime" : job_specs['walltime'],
        "job_mem_per_cpu" : job_specs['mem_per_cpu'], # in MB
        "scaling_function" : job_specs['scaling_fct'],
        "scale_index" : job_specs['scale_index']
        })

where ``path/to/benchmark_dir`` is the path towards the directory where you want your benchmark files to be stored.

Now, at the end of your jobs, a new temporary CSV file will be created in your ``benchmark_dir`` directory, named ``<prog>_<cluster_name>_tmp.csv``, where ``<prog>`` and ``<cluster_name>`` are the names of your program and your cluster, respectively. If the file already exists, a new line will simply be added to it.

Configure the crontab task
--------------------------

Use the ``crontab -e`` command to edit your crontab tasks and add the following line:

.. code-block::

   */15 * * * * bash -l -c "/path/to/cron_benchmark.sh <prog> <cluster_name>" >> path/to/benchmark_dir/crontab_<prog>_<cluster_name>.log 2>&1

where

- 15 is the number of minutes between two consecutive executions of this command (feel free to adjust it at will).
- ``/path/to/cron_benchmark.sh`` is the path towards the crontab script.
- ``crontab_<prog>_<cluster_name>.log`` is a log file that will contain the output of the execution of this crontab script.

Don't forget to also make the ``cron_benchmark.sh`` script executable by entering the following command in your terminal:

.. code-block:: shell

   chmod u+x /path/to/cron_benchmark.sh

Configure the crontab script
----------------------------

In the ``cron_benchmark.sh`` script itself, at the beginning of the file, you will need to specify the path to your benchmark directory:

.. code-block:: shell

   benchmark_dir="path/to/benchmark_dir"

and you will need to **load your Python distribution** (if it is not loaded by default in your user profile configuration).

Now, every 15 minutes, the ``cron_benchmark.sh`` script will be executed. If there is a file named ``<prog>_<cluster_name>_tmp.csv`` in your benchmark directory, it wil archive it into an ``archive`` subdirectory and rename it with the current date. It will then execute ``benchmark.py`` on that file. 

That last script will either create or update the final csv file, named ``<prog>_<cluster_name>.csv`` and placed inside the benchmark directory. The log file of this Python execution can be found in a ``bench_logs`` subdirectory, named ``<prog>_<cluster_name>_<current_date>.log``.
