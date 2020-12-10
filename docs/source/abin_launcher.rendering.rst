***********************
Rendering the templates
***********************

What is rendering?
==================

Rendering the templates consists in creating the input files and job script for our calculations, based on their Jinja templates. The Jinja language offers a relatively easy and intuitive way of creating those templates and allows for more flexibility and adaptability in ``ABIN LAUNCHER``.

The rendering process makes use of three main elements:

- The Jinja templates, defining how the information is presented in the input files and the job script
- The content needed to "fill" those templates with the information specific to each calculation, partially provided by a new YAML file, called the configuration file
- The rendering function that will make the link between those two, i.e. filling the templates with the content

Those elements are presented in more details in the following subsections. If you haven't consulted it yet, it is strongly suggested to take a look at the :doc:`basics.know` section, which contains tutorials that might prove helpful in your understanding of this section.

.. Important::

   This part of the process must almost entirely be **defined by the user themselves**: you will have to create your templates, your configuration files and your rendering functions in order to deal with your specific problematic. This may look a bit bothersome, but it is also what makes the great flexibility of ``ABIN LAUNCHER``, and you will only need to set it up once.

Jinja templates
===============

Rather than trying to explain how to create your templates in a vacuum, let's consider two basic examples to illustrate this process.

.. warning::

   All the Jinja templates must be placed inside the ``templates`` directory of ``ABIN LAUNCHER``.

ORCA input file template
------------------------

First, we will create an input file template for ORCA_. Consider the geometry optimization of |C2H6| through a DFT calculation with the B3LYP functional and the def2-SVP basis set. Based on ORCA's manual_, this is what the input file must be:

.. |C2H6| replace:: C\ :sub:`2`\ H\ :sub:`6`\ 

.. code-block:: text

   ! B3LYP def2-SVP Opt 
   *xyz 0 1
   C         -5.27061        0.53724       -0.00000
   C         -3.75887        0.50632       -0.00000
   H         -5.66181        0.20733       -0.96702
   H         -5.63432        1.55129        0.19093
   H         -5.66857       -0.12330        0.77610
   H         -3.36090        1.16686       -0.77610
   H         -3.39516       -0.50772       -0.19093
   H         -3.36766        0.83624        0.96702
   *

To make a Jinja template for such an input file, we need to replace every variable part of this input by a Jinja variable:

.. literalinclude:: abin_sample/sample_files/sample_orca.inp.jinja
   :language: jinja

We can intuitively tell that:

- ``{{ method }}`` and ``{{ basis_set }}`` replaced *B3LYP* and *def2-SVP*, respectively.
- ``{{ job_type }}`` replaced *Opt*.
- ``{{ charge }}`` and ``{{ multiplicity }}`` replaced the *0* and *1*, respectively.

For the atomic coordinates, since it is impossible to know in advance the number of lines they will contain, we simply use a *for* loop that iterates over each element of the ``coordinates`` list.

That's it! Now we have a basic functioning ORCA input file template. The only thing we will have to do in the future is to pass the content of the variables to Jinja and create the corresponding input file. As long as we know what the input file must look like, we can define a template for pretty much any ab initio program.

Job script template
-------------------

The input file is not the only file that needs to be created. We also need to create the job script, that will give commands to the job scheduler to run the ORCA calculation on the cluster. Let's go back to our example and consider that we are using SLURM as a job scheduler, this is what the job script might look like:

.. code-block:: slurm

   #!/bin/bash

   #SBATCH --output=slurm_output.log
   #SBATCH --job-name=sample_orca_c2h6_svp
   #SBATCH --mail-user=your@email.com
   #SBATCH --mail-type=FAIL
   #SBATCH --time=0-05:00:00
   #SBATCH --ntasks=4
   #SBATCH --mem-per-cpu=500
   #SBATCH --partition=your-cluster-partition

   echo -e "\n================= ORCA execution begins now =================="

   module load ORCA/4.1.0-OpenMPI-3.1.3
   /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca c2h6.inp > c2h6.out 

   echo -e "\n=================  ORCA execution ends now  =================="

Now, we need to replace every variable part of this input by a Jinja variable:

.. code-block:: jinja

   #!/bin/bash

   #SBATCH --output=slurm_output.log
   #SBATCH --job-name={{ profile }}_{{ mol_name }}_{{ config_name }}
   #SBATCH --mail-user={{ user_email }}
   #SBATCH --mail-type={{ mail_type }}
   #SBATCH --time={{ job_walltime }}
   #SBATCH --ntasks={{ job_cores }}
   #SBATCH --mem-per-cpu={{ job_mem_per_cpu }}
   {% if partition != None %}
   #SBATCH --partition={{ partition }}
   {% endif %}

   echo -e "\n================= ORCA execution begins now =================="

   module load ORCA/4.1.0-OpenMPI-3.1.3
   /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca {{ mol_name }}.inp > {{ mol_name }}.out 

   echo -e "\n=================  ORCA execution ends now  =================="

Note that since the partition value is optional, we need to check if it has been specified before putting it in the instructions. Otherwise, we might end up with a ``partition=None`` which won't be recognized by the cluster.

We are almost done, but there is one more thing we need to take into consideration here: ``ABIN LAUNCHER`` puts a copy of the geometry file inside the job subdirectory (see :ref:`abin_out_dir_struct`). If we don't want that file to be overwritten by the new optimized geometry, we need to rename it before ORCA starts running:

.. code-block:: jinja

   #!/bin/bash

   #SBATCH --output=slurm_output.log
   #SBATCH --job-name={{ profile }}_{{ mol_name }}_{{ config_name }}
   #SBATCH --mail-user={{ user_email }}
   #SBATCH --mail-type={{ mail_type }}
   #SBATCH --time={{ job_walltime }}
   #SBATCH --ntasks={{ job_cores }}
   #SBATCH --mem-per-cpu={{ job_mem_per_cpu }}
   {% if partition != None %}
   #SBATCH --partition={{ partition }}
   {% endif %}

   echo -e "Renaming the original .xyz file to avoid overwriting it with the new one."
   cd $SLURM_SUBMIT_DIR
   mv {{ mol_name }}.xyz {{ mol_name }}_ori.xyz

   echo -e "\n================= ORCA execution begins now =================="

   module load ORCA/4.1.0-OpenMPI-3.1.3
   /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca {{ mol_name }}.inp > {{ mol_name }}.out 

   echo -e "\n=================  ORCA execution ends now  =================="

Once again, that's it! Now we have a basic functioning job script template for ORCA jobs. 

However, we can also take this template one step further:

.. literalinclude:: abin_sample/sample_files/sample_orca_job.sh.jinja
   :language: jinja

where the module to load has been replaced by ``{{ set_env }}`` and the command to execute the program has been replaced by ``{{ command }}``. With this, it becomes possible to load and run the same program on another SLURM cluster, where the module and/or the command might be different.

As long as we know what the job script must look like, we can define a template for pretty much any job on any cluster with any job scheduler.

YAML configuration file
=======================

Some of the content for the Jinja variables has already been defined by ``ABIN LAUNCHER`` at this point, with the help of the :doc:`scanning <abin_launcher.scan>` and :doc:`scaling <abin_launcher.job_scale>` processes. For the rest, another YAML file needs to be created, called the configuration file, which contains all the parameters specific to each calculation. Note that just like you can define multiple geometry files, you can define multiple configuration files if the need arises.

To illustrate this part of the process, let's continue with our |C2H6| example. We need the following information:

ORCA input file:

- ``{{ method }}``
- ``{{ basis_set }}``
- ``{{ job_type }}``
- ``{{ charge }}``
- ``{{ multiplicity }}`` 
- ``{{ coordinates }}`` - already defined in the ``file_data`` variable, as built by the :doc:`scanning function <abin_launcher.scan>`

Job script:

- ``{{ profile }}`` - already given as a :ref:`command line argument <abin_arguments>`
- ``{{ mol_name }}`` - name of the geometry file (minus the extension)
- ``{{ config_name }}`` - name of the configuration file (minus the extension)
- ``{{ user_email }}``
- ``{{ mail_type }}``
- ``{{ job_walltime }}`` - already defined by the :doc:`job scaling <abin_launcher.job_scale>` process
- ``{{ job_cores }}`` - already defined by the :doc:`job scaling <abin_launcher.job_scale>` process
- ``{{ job_mem_per_cpu }}`` - already defined by the :doc:`job scaling <abin_launcher.job_scale>` process
- ``{{ partition }}`` - already defined by the :doc:`job scaling <abin_launcher.job_scale>` process
- ``{{ set_env }}``
- ``{{ command }}``

We can then define the configuration file content as:

.. literalinclude:: abin_sample/sample_files/svp.yml
   :language: yaml

For the ``{{ set_env }}`` and ``{{ command }}`` variables, since they are dependent on the cluster, it makes more sense to define them in the :ref:`clusters configuration file <abin_clusters_file>`:

.. code-block:: yaml

   myclusterA:
     profiles:
       sample_orca:
         set_env: module load ORCA/4.1.0-OpenMPI-3.1.3
         command: /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca
   
   myclusterB:
     profiles:
       sample_orca:
         set_env: another-module
         command: another-command

Note that for intuitiveness purposes, the name of the YAML keys is identical to the name of the Jinja variables (``method`` for ``{{ method }}``, ``basis_set`` for ``{{ basis_set }}``, etc.), but nothing prevents you from doing things differently.

.. Caution::

   If for whatever reason, you don't want or need to assign a value to one of the YAML keys in the configuration files (perhaps because the default value is what you need), you must still specify its value as an empty string (``""``). Otherwise, the value will be printed as ``None``, which might cause some problems when running your program.

.. _rendering_fct:

Rendering functions
===================

The Jinja templates define how the information is presented in the input files and the job script. The configuration file, among others, defines the content specific to each calculation. To link the two, ``ABIN LAUNCHER`` calls a function defined in the ``renderer.py`` file, called a **rendering function**, that creates the input files and the job script by filling the templates with their specific content.

General definition
------------------

All the rendering functions must be defined in the ``renderer.py`` file and need to obey some restrictions in order to be callable by ``ABIN LAUNCHER``:

- They take six dictionaries as arguments: ``mendeleev``, ``clusters_cfg``, ``config``, ``file_data``, ``job_specs`` and ``misc`` (see the next subsection for details).
- They must return two variables : 

   - A dictionary where each key is the name of the file and each associated value the rendered content of that file, for example {<input_filename>:<input_content> ; <job_script_name>:<job_script_content>}.
   - The name of the rendered job script, needed by the main script to launch the job on the cluster.

If a problem arises when rendering the templates, an ``AbinError`` exception should be raised with a proper error message (see :ref:`how to handle errors <abin_errors>` for more details).

The six arguments
-----------------

As said in the previous subsection, the rendering functions take six dictionaries as arguments. Since those functions might want various information depending on each specific case, we tried to include as many pertinent details as you might want to refer to during your rendering process. Thus, the six dictionaries are defined as follows:

- ``mendeleev`` is the content of the ``mendeleev.yml`` file.
- ``clusters_cfg`` is the content of the :ref:`clusters configuration file <abin_clusters_file>`.
- ``config`` is the content of the YAML configuration file.
- ``file_data`` is the variable built by the :doc:`scanning function <abin_launcher.scan>`.

- ``job_specs`` contains the information about the resources requirements, defined by the :doc:`job scaling <abin_launcher.job_scale>` process, as well as other details about the job:

   - ``profile``, the name of the profile as it appears in the :ref:`clusters configuration file <abin_clusters_file>` and as it was given :ref:`in the command line <abin_arguments>`.
   - ``scaling_fct``, the name of the chosen :ref:`scaling function <scaling_fcts>`.
   - ``scale_index``, the computed value of the scale_index.
   - ``cluster_name``, the name of the cluster on which ``ABIN LAUNCHER`` is running, as it was given :ref:`in the command line <abin_arguments>`.
   - ``scale_label``, ``scale_limit``, ``partition``, ``walltime``, ``cores`` and ``mem_per_cpu`` are all the information associated to the chosen :ref:`job scale <job_scales>`.

- ``misc`` contains all other pertinent details:

   - ``code_dir`` is the path towards the ``ABIN LAUNCHER`` directory.
   - ``templates_dir`` is the path towards the ``templates`` directory in the ``ABIN LAUNCHER`` directory.
   - ``mol_name`` is the name of the geometry file (minus the extension).
   - ``mol_content`` is the content of the geometry file (in the form of a list where each element is a line of the file).
   - ``config_name`` is the name of the configuration file.

Using the Jinja templates
-------------------------

In order to use and render the Jinja templates, the rendering functions call another function, named ``jinja_render`` and defined in the ``renderer.py`` file. This function is the one that makes the link between the templates and their specific content:

.. autofunction:: renderer.jinja_render

This function receives three key arguments describing where to find the template and what is the corresponding content. It then returns the content of the rendered file in a variable (``output_text``), that will be printed in a file by the main script ``abin_launcher.py``.

Simple function model
---------------------

In essence, the rendering functions have a pretty simple structure: their task is to define the name of each Jinja template and define the needed content for that template (stored in ``<file>_render_vars``). Here is a basic rendering function model:

.. code-block:: python

   def profile_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):

      # Define the names of the templates

      template_input = "name_of_input_template"
      template_script = "name_of_job_script_template"

      # Define the names of the rendered files

      rendered_input = "name_of_created_input_file"
      rendered_script = "name_of_created_job_script_file"

      # Initialize the dictionary that will be returned by the function

      rendered_content = {}

      # Render the template for the input file

      input_render_vars = {
         "jinja_variable1" : value,
         "jinja_variable2" : value,
         "jinja_variable3" : value,
         ...
      }

      rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

      # Render the template for the job script

      script_render_vars = {
         "jinja_variable1" : value,
         "jinja_variable2" : value,
         "jinja_variable3" : value,
         ...
      }

      rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

      # Return the content of the rendered files and the name of the rendered job script

      return rendered_content, rendered_script

A basic example of a rendering function is presented at the end of this section. Note however that some additional steps might be required depending on each specific case.

Calling your rendering function
-------------------------------

The rendering function that will be called by ``ABIN LAUNCHER`` is the one associated with the ``rendering_function`` YAML key defined in the :ref:`clusters configuration file <abin_clusters_file>`:

.. code-block:: yaml

   mycluster:
      profiles:
         myprofile1:
            rendering_function: name-of-rendering-function
         myprofile2:
            rendering_function: name-of-rendering-function

where ``mycluster`` corresponds to the name of your cluster (given as a :ref:`command line argument <abin_arguments>`) while ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to run (such as ``orca`` or ``qchem``). This way, a different rendering function can be assigned to each profile.

Example of a rendering function
-------------------------------

Let's end this section by making the rendering function associated with our geometry optimization example.

Review of the template and configuration files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before defining the function, we need to review our different files:

First, we have the input file template:

.. literalinclude:: abin_sample/sample_files/sample_orca.inp.jinja
   :language: jinja
   :caption: sample_orca.inp.jinja

Then, the job script template:

.. literalinclude:: abin_sample/sample_files/sample_orca_job.sh.jinja
   :language: jinja
   :caption: sample_orca_job.sh.jinja

Third, the configuration file:

.. literalinclude:: abin_sample/sample_files/svp.yml
   :language: yaml
   :caption: svp.yml

Finally, let's consider that we have the following clusters configuration file:

.. literalinclude:: abin_sample/sample_files/clusters.yml
   :language: yaml
   :caption: clusters.yml

where ``lemaitre3`` is the name of our cluster.

Function definition
~~~~~~~~~~~~~~~~~~~

Now that we know what our different files look like, let's define the rendering function, called ``sample_orca_render``, associated with our example:

.. code-block:: python

   def sample_orca_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):

We first have to define the names of our Jinja templates.

.. code-block:: python

    template_input = "sample_orca.inp.jinja"
    template_script = "sample_orca_job.sh.jinja"

Then we have to define the names of our rendered files.

.. code-block:: python

    rendered_input = misc['mol_name'] + ".inp"
    rendered_script = "orca_job.sh"

Next comes the initialization of the ``rendered_content`` dictionary, where each key will be the name of the file and each value the rendered content of that file. (This is what the function will return.)

.. code-block:: python

    rendered_content = {}

Now, let's define the content of our input file. This content is stored in the ``input_render_vars`` dictionary, where each key corresponds to one of our Jinja variables. This dictionary is the key part of the rendering function as it defines the link between the Jinja template and the various information needed to fill it.

.. code-block:: python

      input_render_vars = {
         "method" : config['method'],
         "basis_set" : config['basis_set'],
         "job_type" : config['job_type'],
         "charge" : config['charge'],
         "multiplicity" : config['multiplicity'],
         "coordinates" : file_data['atomic_coordinates']
      }

Our first template is now ready to be rendered, let's call the ``jinja_render`` function and store the result into ``rendered_content`` under the key ``rendered_input`` (the name of our rendered file):

.. code-block:: python

    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, input_render_vars)

Let's proceed with the second template in the same manner:
     
.. code-block:: python

      script_render_vars = {  
         "mol_name" : misc['mol_name'],
         "config_name" : misc['config_name'],
         "user_email" : config['user_email'],
         "mail_type" : config['mail_type'],
         "job_walltime" : job_specs['walltime'],
         "job_cores" : job_specs['cores'],
         "job_mem_per_cpu" : job_specs['mem_per_cpu'],
         "partition" : job_specs['partition'],     
         "set_env" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['set_env'],       
         "command" : clusters_cfg[job_specs['cluster_name']]['profiles'][job_specs['profile']]['command'],
         "profile" : job_specs['profile']
      }

      rendered_content[rendered_script] = jinja_render(misc['templates_dir'], template_script, script_render_vars)

Finally, we just need to return ``rendered_content`` and ``rendered_script`` to the main script, so that the rendered content can be printed to the different files and the job can be launched on the cluster:

.. code-block:: python

    return rendered_content, rendered_script

Our function is now ready. This is what it ends up looking like with proper comments and documentation:

.. literalinclude:: abin_sample/sample_files/renderer.py
   :language: python
   :lines: 46-132

.. Hyperlink targets

.. _manual: https://sites.google.com/site/orcainputlibrary/generalinput
.. _ORCA: https://www.faccts.de/orca/