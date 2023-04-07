***********************
Rendering the templates
***********************

What is rendering?
==================

Rendering the templates consists in creating the parameters file and job script for our QOCT-GRAD calculations, based on their Jinja templates. The Jinja language offers a relatively easy and intuitive way of creating those templates and allows for more flexibility and adaptability in ``CONTROL LAUNCHER``.

The rendering process makes use of three main elements:

- The Jinja templates, defining how the information is presented in the parameters file and the job script
- The content needed to "fill" those templates with the information specific to each calculation, partially provided by a new YAML file, called the configuration file
- The rendering function that will make the link between those two, i.e. filling the templates with the content

Those elements are presented in more details in the following subsections. If you haven't consulted it yet, it is strongly suggested to take a look at the :doc:`basics.know` section, which contains tutorials that might prove helpful in your understanding of this section.

Jinja templates
===============

Rather than trying to explain how to create your templates in a vacuum, let's consider two basic examples to illustrate this process.

.. warning::

   All the Jinja templates must be placed inside the ``templates`` directory of ``CONTROL LAUNCHER``.

Parameters file template
------------------------

.. todo::

   COMING SOON

Job script template
-------------------

.. todo::

   COMING SOON

YAML configuration file
=======================

.. todo::

   COMING SOON

.. _control_rendering_fct:

Rendering functions
===================

The Jinja templates define how the information is presented in the parameters file and the job script. The configuration file, among others, defines the content specific to each calculation. To link the two, ``CONTROL LAUNCHER`` calls a function defined in the ``control_renderer.py`` file, called a **rendering function**, that creates the parameters file and the job script by filling the templates with their specific content.

.. Important::

   The rendering function will be called for each transition-configuration combination, i.e. for each transition defined by the :doc:`transition function <control_launcher.transitions>` coupled with each configuration file.

General definition
------------------

All the rendering functions must be defined in the ``control_renderer.py`` file and need to obey some restrictions in order to be callable by ``CONTROL LAUNCHER``:

- They take six dictionaries as arguments: ``clusters_cfg``, ``config``, ``system``, ``data``, ``job_specs`` and ``misc`` (see the next subsection for details).
- They must return two variables : 

   - A dictionary where each key is the name of the file and each associated value the rendered content of that file, for example {<parameters_filename>:<parameters_content> ; <job_script_name>:<job_script_content>}.
   - The name of the rendered job script, needed by the main script to launch the job on the cluster.

If a problem arises when rendering the templates, a ``ControlError`` exception should be raised with a proper error message (see :ref:`how to handle errors <control_errors>` for more details).

The six arguments
-----------------

As said in the previous subsection, the rendering functions take six dictionaries as arguments. Since those functions might want various information depending on each specific case, we tried to include as many pertinent details as you might want to refer to during your rendering process. Thus, the six dictionaries are defined as follows:

- ``clusters_cfg`` is the content of the :ref:`clusters configuration file <control_clusters_file>`.
- ``config`` is the content of the YAML configuration file.
- ``system`` is the variable built by the :doc:`system modelling <control_launcher.modelling>` process.
- ``data`` contains the paths towards the ``data`` directory and its content, created during the :doc:`system modelling <control_launcher.modelling>` process and when :doc:`determining the transitions <control_launcher.transitions>`:

   - ``main_path`` is the absolute path towards the the ``data`` directory, created during the :doc:`system modelling <control_launcher.modelling>` process.
   - ``mime_path``, ``momdip_mtx_path``, ``eigenvalues_path``, ``eigenvectors_path``, ``transpose_path`` and ``momdip_es_mtx_path`` are the path towards each of the file created during the :doc:`system modelling <control_launcher.modelling>` process.
   - ``init_path`` and ``target_path`` are the path towards the initial and target states file of the considered transition (defined by the :doc:`transition function <control_launcher.transitions>` and placed inside the ``data`` directory).

- ``job_specs`` contains the information about the resources requirements, defined by the :doc:`job scaling <control_launcher.job_scale>` process, as well as other details about the job:

   - ``profile`` is the name of the profile as it appears in the :ref:`clusters configuration file <control_clusters_file>` and as it was given :ref:`in the command line <control_arguments>`.
   - ``scale_index`` is the computed value of the scale_index.
   - ``cluster_name`` is the name of the cluster on which ``CONTROL LAUNCHER`` is running, as it was given :ref:`in the command line <control_arguments>`.
   - ``scale_label``, ``scale_limit``, ``partition``, ``walltime`` and ``memory`` are all the information associated to the chosen :ref:`job scale <control_job_scales>`.

- ``misc`` contains all other pertinent details:

   - ``code_dir`` is the path towards the ``CONTROL LAUNCHER`` directory.
   - ``templates_dir`` is the path towards the ``templates`` directory in the ``CONTROL LAUNCHER`` directory.
   - ``source_name`` is the name of the source file (minus a possible extension).
   - ``source_content`` is the content of the source file (in the form of a list where each element is a line of the file).
   - ``mol_dir`` is the path towards the molecule directory, bearing the name of the source file and containing the ``data`` directory and the job subdirectories.
   - ``config_name`` is the name of the configuration file.
   - ``transition_label`` is the label of the considered transition.
   - ``transitions_list``, the list of dictionaries built when :doc:`determining the transitions <control_launcher.transitions>`.

Using the Jinja templates
-------------------------

In order to use and render the Jinja templates, the rendering functions call another function, named ``jinja_render`` and defined in the ``control_renderer.py`` file. This function is the one that makes the link between the templates and their specific content:

.. autofunction:: control_renderer.jinja_render

This function receives three key arguments describing where to find the template and what is the corresponding content. It then returns the content of the rendered file in a variable (``output_text``), that will be printed in a file by the main script ``control_launcher.py``.

Simple function model
---------------------

In essence, the rendering functions have a pretty simple structure: their task is to define the name of each Jinja template and define the needed content for that template (stored in ``<file>_render_vars``). Here is a basic rendering function model:

.. code-block:: python

   def my_rendering_function(clusters_cfg:dict, config:dict, system:dict, data:dict, job_specs:dict, misc:dict):

      # Define the names of the templates

      template_param = "name_of_parameters_template"
      template_script = "name_of_job_script_template"

      # Define the names of the rendered files

      rendered_param = "name_of_created_parameters_file"
      rendered_script = "name_of_created_job_script_file"

      # Initialize the dictionary that will be returned by the function

      rendered_content = {}

      # Render the template for the parameters file

      param_render_vars = {
         "jinja_variable1" : value,
         "jinja_variable2" : value,
         "jinja_variable3" : value,
         ...
      }

      rendered_content[rendered_param] = jinja_render(misc['templates_dir'], template_param, param_render_vars)

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

The rendering function that will be called by ``CONTROL LAUNCHER`` is the one associated with the ``rendering_function`` YAML key defined in the :ref:`clusters configuration file <control_clusters_file>`:

.. code-block:: yaml

   mycluster:
      profiles:
         myprofile1:
            rendering_function: name-of-rendering-function
         myprofile2:
            rendering_function: name-of-rendering-function

where ``mycluster`` corresponds to the name of your cluster (given as a :ref:`command line argument <control_arguments>`) while ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to use. This way, a different rendering function can be assigned to each profile.

Example of a rendering function
-------------------------------

.. todo::

   COMING SOON

Review of the template and configuration files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before defining the function, we need to review our different files:

First, we have the parameters file template:

.. literalinclude:: control_sample/sample_param.nml.jinja
   :language: jinja
   :caption: sample_param.nml.jinja

Then, the job script template:

.. literalinclude:: control_sample/sample_qoctra_job.sh.jinja
   :language: jinja
   :caption: sample_qoctra_job.sh.jinja

Third, the configuration file:

.. literalinclude:: control_sample/config.yml
   :language: yaml
   :caption: config.yml

Finally, let's consider that we have the following clusters configuration file:

.. literalinclude:: control_sample/clusters.yml
   :language: yaml
   :caption: clusters.yml

where ``lemaitre3`` is the name of our cluster.

Function definition
~~~~~~~~~~~~~~~~~~~

.. todo::

   COMING SOON

.. Hyperlink targets
