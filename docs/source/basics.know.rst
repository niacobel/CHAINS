*********************
What You Need To Know
*********************

About YAML
==========

YAML_ is a "human friendly data serialization standard for all programming languages". (quote from the official website)

Here is a small youtube tutorial made by Mike Dane that will cover the basics of the YAML syntax:

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/cdLNKUoMc6c" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

.. .. youtube:: cdLNKUoMc6c

|
| With that, you know everything you need to know about YAML in order to use CHAINS! YAML is used here for the configuration files, containing all the parameters needed for the different scripts and calculations. As you may have seen in the video, the way those files are structured makes them very easy to use and to understand.

About Jinja2
============

"Jinja_ is a modern and designer-friendly templating language for Python, modelled after Django’s templates. It is fast, widely used and secure with the optional sandboxed template execution environment." (quote from the official website)

Here is a small youtube tutorial made by Jason Rigden that will cover the basics of the Jinja2 templating engine:

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/bxhXQG1qJPM" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

.. .. youtube:: bxhXQG1qJPM

|
| With that, you know everything you need to know about Jinja2 in order to use CHAINS! We use the Jinja language to create templates for the input files needed by ORCA_, Q-CHEM_ and QOCT-RA_. We also use it for the job scripts, containing all instructions for the job scheduler (SLURM_ for the CECI clusters). As you may have seen in the video, the Jinja language offers a relatively easy and intuitive way of creating those templates, which makes it easier to adapt the templates to pretty much any input files for any ab initio program such as Gaussian_ or Molpro_, and any other job scheduler such as Torque_.

.. todo::
   COMING SOON
   Theoretical background
   Python knowledge
   Clusters and job scheduler

.. Hyperlink targets

.. _Gaussian: https://gaussian.com/
.. _Jinja: https://jinja.palletsprojects.com/en/2.11.x/ 
.. _Molpro: https://www.molpro.net/
.. _ORCA: https://www.faccts.de/orca/
.. _Q-CHEM: https://www.q-chem.com/
.. _QOCT-RA: https://gitlab.com/dynaq.cqp/QOCT-RA
.. _SLURM: https://slurm.schedmd.com/documentation.html
.. _Torque: https://github.com/adaptivecomputing/torque
.. _YAML: https://yaml.org/