***********************
Rendering the templates
***********************

What is rendering?
==================

Rendering the templates consists in creating the input files and job instructions file for our calculations, based on their Jinja templates. The Jinja language offers a relatively easy and intuitive way of creating those templates and allows for more flexibility and adaptability in ``ABIN LAUNCHER``.

The rendering process makes use of three main elements:

- The Jinja templates, defining how the information is presented in the input files and the job instructions file
- The content needed to "fill" those templates with the information specific to each calculation, partially provided by a new YAML file, called the configuration file
- The rendering function that will make the link between those two, i.e. filling the templates with the content

Those elements are presented in more details in the following subsections. 

.. note::

   This part of the process must almost entirely be defined by the user themselves: you will have to create your templates, your configuration files and your rendering functions in order to deal with your specific problematic.

Jinja templates
===============

To illustrate how to create those Jinja templates, let's consider two basic examples. 

.. warning::

   All the Jinja templates must be placed inside the ``Templates`` directory of ``ABIN LAUNCHER``.

ORCA input file template
------------------------

First, we will create an input file template for ORCA_. Consider the geometry optimization of |C2H6| through a DFT calculation with the B3LYP functional and the def2-SVP basis set. Based on ORCA's manual_, this is what the input file must be:

.. |C2H6| replace:: C\ :sub:`2`\ H\ :sub:`6`\ 

.. code-block::

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

.. code-block:: jinja

   ! {{ method }} {{ basis_set }} {{ job_type }} 
   *xyz {{ charge }} {{ multiplicity }}
   {% for coordinate_line in coordinates -%}
   {{ coordinate_line }}
   {% endfor -%}
   *

We can intuitively tell that:

- ``{{ method }}`` and ``{{ basis_set }}`` replaced *B3LYP* and *def2-SVP*, respectively.
- ``{{ job_type }}`` replaced *Opt*.
- ``{{ charge }}`` and ``{{ multiplicity }}`` replaced the *0* and *1*, respectively.

For the atomic coordinates, since it is impossible to know in advance the number of lines they will contain, we simply use a *for* loop that iterates over each element of the ``coordinates`` list.

That's it! Now we have a basic functioning ORCA input file template. The only thing we will have to do in the future is to pass the content of the variables to Jinja and create the corresponding input file. As long as we know what the input file must look like, we can define a template for pretty much any ab initio program.

Job instructions file template
------------------------------

The input file is not the only file that needs to be created. We also need to create the job instructions file, that will give commands to the job scheduler to run the ORCA calculation on the cluster. Let's go back to our example and consider that we are using SLURM as a job scheduler, this is what the job instructions file might look like:

.. code-block:: shell

   #!/bin/bash

   #SBATCH --output=slurm_output.log
   #SBATCH --job-name=c2h6_orca
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
   #SBATCH --job-name={{ mol_name }}_orca
   #SBATCH --mail-user={{ user_email }}
   #SBATCH --mail-type={{ mail_type }}
   #SBATCH --time={{ job_walltime }}
   #SBATCH --ntasks={{ job_cores }}
   #SBATCH --mem-per-cpu={{ job_mem_per_cpu }}
   #SBATCH --partition={{ partition }}

   echo -e "\n================= ORCA execution begins now =================="

   module load ORCA/4.1.0-OpenMPI-3.1.3
   /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca {{ mol_name }}.inp > {{ mol_name }}.out 

   echo -e "\n=================  ORCA execution ends now  =================="

Once again, that's it! Now we have a basic functioning job instructions file template for ORCA jobs. However, we can also take this template one step further:

.. code-block:: jinja

   #!/bin/bash

   #SBATCH --output=slurm_output.log
   #SBATCH --job-name={{ mol_name }}_{{ prog }}
   #SBATCH --mail-user={{ user_email }}
   #SBATCH --mail-type={{ mail_type }}
   #SBATCH --time={{ job_walltime }}
   #SBATCH --ntasks={{ job_cores }}
   #SBATCH --mem-per-cpu={{ job_mem_per_cpu }}
   #SBATCH --partition={{ partition }}

   echo -e "\n================= {{ prog }} execution begins now =================="

   {{ set_env }}
   {{ command }} {{ mol_name }}.inp > {{ mol_name }}.out 

   echo -e "\n=================  {{ prog }} execution ends now  =================="

where the module to load has been replaced by ``{{ set_env }}`` and the command to execute the program has been replaced by ``{{ command }}``.

With this, it becomes possible to

- load and run another program.
- load and run the same program on another SLURM cluster, where the module and/or the command might be different.

As long as we know what the job instructions file must look like, we can define a template for pretty much any job.

YAML configuration file
=======================

Some of the content for the Jinja variables has already been defined by ``ABIN LAUNCHER`` at this point, with the help of the :doc:`scanning <abin_launcher.scan>` and :doc:`scaling <abin_launcher.job_scale>` processes. For the rest, another YAML file needs to be created, called the configuration file, which contains all the parameters specific to each calculation. Note that just like you can define multiple geometry files, you can define multiple configuration file if the need arises.

To illustrate this part of the process, let's continue with our |C2H6| example. We need the following information:

ORCA input file:

- ``{{ method }}``
- ``{{ basis_set }}``
- ``{{ job_type }}``
- ``{{ charge }}``
- ``{{ multiplicity }}`` 
- ``{{ coordinates }}`` - already defined in the ``file_data`` variable, as built by the :doc:`scanning function <abin_launcher.scan>`

Job instructions file:

- ``{{ mol_name }}`` - defined as the same name as the geometry file (minus the extension)
- ``{{ prog }}`` - already given as a :ref:`command line argument <abin_arguments>`
- ``{{ user_email }}``
- ``{{ mail_type }}``
- ``{{ job_walltime }}`` - already defined by the :doc:`job scaling <abin_launcher.job_scale>` process
- ``{{ job_cores }}`` - already defined by the :doc:`job scaling <abin_launcher.job_scale>` process
- ``{{ job_mem_per_cpu }}`` - already defined by the :doc:`job scaling <abin_launcher.job_scale>` process
- ``{{ partition }}`` - already defined by the :doc:`job scaling <abin_launcher.job_scale>` process
- ``{{ set_env }}``
- ``{{ command }}``

We can then define the configuration file content as:

.. code-block:: yaml

   method: B3LYP
   basis-set: def2-SVP
   job-type: Opt
   charge: 0
   multiplicity: 1
   user-email: your@email.com
   mail-type: FAIL

For the ``{{ set_env }}`` and ``{{ command }}`` variables, since they are dependent on the cluster, it makes more sense to define them in the :ref:`clusters configuration file <clusters_file>`:

.. code-block:: yaml

   myclusterA:
     progs:
       orca:
         set_env: module load ORCA/4.1.0-OpenMPI-3.1.3
         command: /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca
   
   myclusterB:
     progs:
       orca:
         set_env: another-module
         command: another-command

Note that for intuitivity purposes, the name of the YAML keys is close to identical to the name of the Jinja variables (``method`` for ``{{ method }}``, ``basis-set`` for ``{{ basis_set }}``, etc.), but nothing prevents you from doing things differently.

.. _rendering_fct:

Rendering functions
===================

The Jinja templates define how the information is presented in the input files and the job instructions file. The configuration file, among others, defines the content specific to each calculation. To link the two, ``ABIN LAUNCHER`` calls a function defined in the ``renderer.py`` file, called a **rendering function**, that creates the input files and the job instructions file by filling the templates with the specific content.

Using the Jinja templates
-------------------------

The ``jinja_render`` function in the ``renderer.py`` file defines how to interact with the Jinja templates:

.. autofunction:: renderer.jinja_render

This function receives three key arguments describing where to find the template and what is the corresponding content. It then returns the content of the rendered file in a variable (``output_text``), that will be printed in a file by the main script ``abin_launcher.py``.

General definition
------------------

All the rendering functions must be defined in the ``renderer.py`` file and need to obey some restrictions in order to be callable by ``ABIN LAUNCHER``:

- They need to be called *prog_render*, where *prog* is the name of the program as it appears in the :ref:`clusters configuration file <clusters_file>` and as it was given :ref:`in the command line <abin_arguments>`. 
- They take six dictionaries as arguments: ``mendeleev``, ``clusters_cfg``, ``config``, ``file_data``, ``job_specs`` and ``misc`` (see the next subsection for details).
- They must return a dictionary where each key is the name of the file and each value the rendered content of that file, e.g. {name_of_input_file:content_of_input_file ; name_of_job_file:content_of_job_file}.
- The name of the rendered job instructions file must be taken from the ``job_instructions`` key in the :ref:`clusters configuration file <clusters_file>`.

If a problem arises when rendering the templates, an ``AbinError`` exception should be raised with a proper error message (see :ref:`how to handle errors <abin_errors>` for more details).

The six arguments
-----------------

As said in the previous section, the rendering functions take six dictionaries as arguments. Since those functions might want various information depending on each specific case, we tried to include as many pertinent details that you might want to refer to during your rendering process. Thus, the six dictionaries are defined as follows:

- ``mendeleev`` is the content of the ``mendeleev.yml`` file
- ``clusters_cfg`` is the content of the :ref:`clusters configuration file <clusters_file>`
- ``config`` is the content of the YAML configuration file
- ``file_data`` is the variable built by the :doc:`scanning function <abin_launcher.scan>`.

- ``job_specs`` contains the information about the resources requirements, defined by the :doc:`job scaling <abin_launcher.job_scale>` process, as well as other details about the job:

   - ``prog``, the name of the program as it appears in the :ref:`clusters configuration file <clusters_file>` and as it was given :ref:`in the command line <abin_arguments>`. You can either use this value or explicitly state it in the code since you already know it.
   - ``scaling_fct``, the name of the chosen :ref:`scaling function <scaling_fcts>`
   - ``scale_index``, the computed value of the scale_index
   - ``cluster_name``, the name of the cluster on which ``ABIN LAUNCHER`` is running, as it was given :ref:`in the command line <abin_arguments>`
   - ``scale_label``, ``scale_limit``, ``partition``, ``walltime``, ``cores`` and ``mem_per_cpu`` are all the information associated to the chosen :ref:`job scale <job_scales>`.

- ``misc`` contains all other pertinent details:

   - ``code_dir`` is the path towards the ``ABIN LAUNCHER`` directory
   - ``path_tpl_dir`` is the path towards the ``Templates`` directory in the ``ABIN LAUNCHER`` directory.
   - ``mol_name`` is the name of the geometry file (minus the extension)
   - ``config_name`` is the name of the configuration file

Simple function model
---------------------

In essence, the rendering functions have a pretty simple structure: their task is to define the name of each Jinja template and define the needed content for that template (stored in ``render_vars``). Here is a basic rendering function model:

.. code-block:: python

   def prog_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):

      # Define the names of the templates

      tpl_inp = "name_of_input_template"
      tpl_inst = "name_of_job_instructions_template"

      # Define the names of the rendered files

      rnd_input = "name_of_created_input_file"
      rnd_inst = clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['job_instructions']

      # Initialize the dictionary that will be returned by the function

      rendered_content = {}

      # Render the template for the input file

      render_vars = {
         "jinja_variable1" : value,
         "jinja_variable2" : value,
         "jinja_variable3" : value,
         ...
      }

      rendered_content[rnd_input] = jinja_render(misc['path_tpl_dir'], tpl_inp, render_vars)

      # Render the template for the job instructions file

      render_vars = {
         "jinja_variable1" : value,
         "jinja_variable2" : value,
         "jinja_variable3" : value,
         ...
      }

      rendered_content[rnd_inst] = jinja_render(misc['path_tpl_dir'], tpl_inst, render_vars)

      # Return the content of the rendered files

      return rendered_content

A concrete example of rendering function is presented in the next subsection. Note however that some additional steps might be required depending on each specific case.

Example of a rendering function
-------------------------------

Let's end this section by making the rendering function associated with our geometry optimization example.

Review of the template and configuration files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before defining the function, we need to review our different files:

First, we have the input file template, that will be named ``orca.inp.jinja``:

.. code-block:: jinja

   ! {{ method }} {{ basis_set }} {{ job_type }} 
   *xyz {{ charge }} {{ multiplicity }}
   {% for coordinate_line in coordinates -%}
   {{ coordinate_line }}
   {% endfor -%}
   *

Then, the job instructions file template, that will be named ``orca_job.sh.jinja``:

.. code-block:: jinja

   #!/bin/bash

   #SBATCH --output=slurm_output.log
   #SBATCH --job-name={{ mol_name }}_{{ prog }}
   #SBATCH --mail-user={{ user_email }}
   #SBATCH --mail-type={{ mail_type }}
   #SBATCH --time={{ job_walltime }}
   #SBATCH --ntasks={{ job_cores }}
   #SBATCH --mem-per-cpu={{ job_mem_per_cpu }}
   #SBATCH --partition={{ partition }}

   echo -e "\n================= {{ prog }} execution begins now =================="

   {{ set_env }}
   {{ command }} {{ mol_name }}.inp > {{ mol_name }}.out 

   echo -e "\n=================  {{ prog }} execution ends now  =================="

Third, the configuration file:

.. code-block:: yaml

   method: B3LYP
   basis-set: def2-SVP
   job-type: Opt
   charge: 0
   multiplicity: 1
   user-email: your@email.com
   mail-type: FAIL

Finally, let's consider that we have the following clusters configuration file:

.. code-block:: yaml

   myclusterA: 
     progs:
       orca:
         job_instructions: orca_job.sh
         set_env: module load ORCA/4.1.0-OpenMPI-3.1.3
         command: /opt/cecisw/arch/easybuild/2018b/software/ORCA/4.1.0-OpenMPI-3.1.3/orca           
         scaling_function: total_nb_elec
         job_scales:
           - label: tiny
             scale_limit: 50
             partition_name: default
             time: 0-00:10:00
             cores: 4 
             mem_per_cpu: 500 # in MB
             delay_command:
           - label: small
             scale_limit: 1000
             partition_name: Def
             time: 2-00:00:00
             cores: 4
             mem_per_cpu: 1000 # in MB
             delay_command:
           - label: medium
             scale_limit: 1500
             partition_name: Def
             time: 5-00:00:00
             cores: 8
             mem_per_cpu: 2000 # in MB
             delay_command: --begin=now+60
           - label: big
             scale_limit: 2000
             partition_name: Long
             time: 10-00:00:00
             cores: 16
             mem_per_cpu: 4000 # in MB  
             delay_command: --begin=now+120

Function definition
~~~~~~~~~~~~~~~~~~~

Now that we know what our different files look like, let's define the rendering function, called ``orca_render``, associated with our example:

.. code-block:: python

   def orca_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):

We first have to define the names of our Jinja templates.

.. code-block:: python

    tpl_inp = "orca.inp.jinja"
    tpl_inst = "orca_job.sh.jinja"

Then we have to define the names of our rendered files.

.. code-block:: python

    rnd_input = misc['mol_name'] + ".inp"
    rnd_inst = clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['job_instructions']

Next comes the initialization of the ``rendered_content`` dictionary, where each key will be the name of the file and each value the rendered content of that file. (This is what the function will return.)

.. code-block:: python

    rendered_content = {}

Now, let's define the content of our input file. This content is stored in the ``render_vars`` dictionary, where each key corresponds to one of our Jinja variables. This dictionary is the key part of the rendering function as it makes a direct link between the Jinja templates and the various information needed to fill it.

.. code-block:: python

      render_vars = {
         "method" : config['method'],
         "basis_set" : config['basis-set'],
         "job_type" : config['job-type'],
         "charge" : config['charge'],
         "multiplicity" : config['multiplicity'],
         "coordinates" : file_data['atomic_coordinates']
      }

Our first template is now ready to be rendered, let's call the ``jinja_render`` function and store the result into ``rendered_content`` under the key ``rnd_input`` (the name of our rendered file):

.. code-block:: python

    rendered_content[rnd_input] = jinja_render(misc['path_tpl_dir'], tpl_inp, render_vars)

Let's proceed with the second template in the same manner:
     
.. code-block:: python

      render_vars = {  
         "mol_name" : misc['mol_name'],
         "user_email" : config['user-email'],
         "mail_type" : config['mail-type'],
         "job_walltime" : job_specs['walltime'],
         "job_cores" : job_specs['cores'],
         "job_mem_per_cpu" : job_specs['mem_per_cpu'],
         "partition" : job_specs['partition'],     
         "set_env" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['set_env'],       
         "command" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['command'],
         "prog" : job_specs['prog']
      }

      rendered_content[rnd_inst] = jinja_render(misc['path_tpl_dir'], tpl_inst, render_vars)

Finally, we just need to return ``rendered_content`` to the main script so that the rendered content can be printed to the respective files:

.. code-block:: python

    return rendered_content

Our function is now ready. This is what it ends up looking like with proper comments and documentation:

.. code-block:: python

   def orca_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
       """Renders the job instructions file and the input file associated with the ORCA program.

       Parameters
       ----------
       mendeleev : dict
           Content of AlexGustafsson's Mendeleev Table YAML file (found at https://github.com/AlexGustafsson/molecular-data).
           Unused in this function.
       clusters_cfg : dict
           Content of the YAML clusters configuration file.
       config : dict
           Content of the YAML configuration file.
       file_data : dict
           Information extracted by the scanning function from the geometry file.
       job_specs : dict
           Contains all information related to the job.
       misc : dict
           Contains all the additional variables that did not pertain to the other arguments.
   
       Returns
       -------
       rendered_content : dict
           Dictionary containing the text of all the rendered files in the form of <filename>: <rendered_content>.
       
       Notes
       -----
       Pay a particular attention to the render_vars dictionary, it contains all the definitions of the variables appearing in your jinja template.
       """

      # Define the names of the templates

      tpl_inp = "orca.inp.jinja"
      tpl_inst = "orca_job.sh.jinja"

      # Define the names of the rendered files

      rnd_input = misc['mol_name'] + ".inp"
      rnd_inst = clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['job_instructions']

      # Initialize the dictionary that will be returned by the function

      rendered_content = {}

      # Render the template for the input file

      render_vars = {
         "method" : config['method'],
         "basis_set" : config['basis-set'],
         "job_type" : config['job-type'],
         "charge" : config['charge'],
         "multiplicity" : config['multiplicity'],
         "coordinates" : file_data['atomic_coordinates']
      }

      rendered_content[rnd_input] = jinja_render(misc['path_tpl_dir'], tpl_inp, render_vars)

      # Render the template for the job instructions file

      render_vars = {  
         "mol_name" : misc['mol_name'],
         "user_email" : config['user-email'],
         "mail_type" : config['mail-type'],
         "job_walltime" : job_specs['walltime'],
         "job_cores" : job_specs['cores'],
         "job_mem_per_cpu" : job_specs['mem_per_cpu'],
         "partition" : job_specs['partition'],     
         "set_env" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['set_env'],       
         "command" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['command'],
         "prog" : job_specs['prog']
      }

      rendered_content[rnd_inst] = jinja_render(misc['path_tpl_dir'], tpl_inst, render_vars)

      # Return the content of the rendered files

      return rendered_content

.. Hyperlink targets

.. _manual: https://sites.google.com/site/orcainputlibrary/generalinput
.. _ORCA: https://www.faccts.de/orca/