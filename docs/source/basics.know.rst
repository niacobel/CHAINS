*********************
What you need to know
*********************

In this section, we will cover what programming knowledge you will need in order to use CHAINS or some of its main scripts by introducing you to some great tutorials that can teach you all about it.

About your jobs and your job scheduler
======================================

First things first, it is important to clarify what CHAINS *does not do*. CHAINS cannot help you in figuring out how to run calculations with your programs, or how to build input files for them. CHAINS (and its main scripts) can help you automate that process, but it cannot guess what it is supposed to do. Start by defining what your methodology is and what it implies. Read the documentation of your programs as well as the one for your job scheduler or cluster. Once you know what to do and how to do it manually, it is time to automate it in order to process more jobs more quickly!

Another thing to note is that, since CHAINS was developed on SLURM clusters, our examples will involve SLURM commands. If you use a different job scheduler, you can use the `rosetta PDF <https://slurm.schedmd.com/rosetta.pdf>`_ from the SLURM documentation to "translate" our examples to your case.

About Python
============

Python_ 3.5 is the main language of CHAINS. If you just want to use the scripts without customizing them, a beginner's knowledge of Python should be more than enough. Otherwise, you might need to delve just a bit deeper into that language, although we tried to stick with the basics as much as possible. If some lines of code use anything a bit unusual, a link is usually given in the comments for reference.

In order to learn Python, here is a YouTube tutorials playlist made by Corey Schafer:

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/videoseries?list=PL-osiE80TeTskrapNbzXhwoFUiLCjGgY7" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

| 
| As far as written tutorials go, the Python Guru `website <https://thepythonguru.com/>`_ offers some great ones, as does `TutorialsTeacher <https://www.tutorialsteacher.com/python>`_. The `official documentation <https://docs.python.org/3/tutorial/index.html>`_ is also quite complete. For French speakers, OpenClassRooms constitutes `a great option <https://openclassrooms.com/fr/courses/235344-apprenez-a-programmer-en-python>`_.

Pay a special attention to **Python dictionaries**, as they will be the main focus of the work you will need to do in order to use CHAINS and its main scripts.

About YAML
==========

YAML_ is a "human friendly data serialization standard for all programming languages" (quote from the official website). YAML is used here for the configuration files, containing all the parameters needed for the different scripts and calculations.

Here is a small YouTube tutorial made by Mike Dane that will cover the basics of the YAML syntax:

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/cdLNKUoMc6c" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

.. .. youtube:: cdLNKUoMc6c

|
| When reading the YAML files, PyYAML_ transforms them into Python dictionaries. For example, let's consider the following YAML file:

.. code-block:: yaml
   :caption: example.yml

   Name: Bobby Smith
   Age: 32
   Hobbies:
     - Basketball
     - Driving
     - Watching movies
   Family:
     Wife: Johanna
     Son: John

and the Python script that will read it:

.. code-block:: python
   :caption: read.py

   import yaml

   with open('example.yml') as f:
      data = yaml.load(f, Loader=yaml.FullLoader)
   
   print(data)

When executing the python script, we can see that the ``data`` variable is a Python dictionary:

.. code-block:: shell

   $ python read.py
   {'Name': 'Bobby Smith', 'Age': 32, 'Hobbies': ['Basketball', 'Driving', 'Watching movies'], 'Family': {'Wife': 'Johanna', 'Son': 'John'}}

As you can see, the ``Hobbies`` values have become a Python list and the nested values in ``Family`` have become a Python dictionary. If we wanted to access the name of Bobby's wife, we could use:

.. code-block:: python

   wife_name = data['Family']['Wife']

With that, you know everything you need to know about YAML in order to use CHAINS! As you may have seen in the video, the way those files are structured makes them very easy to use and to understand.

About Jinja2
============

"Jinja_ is a modern and designer-friendly templating language for Python, modelled after Django’s templates. It is fast, widely used and secure with the optional sandboxed template execution environment" (quote from the official website). We use the Jinja language to create templates for the input files needed by Gaussian_, Q-CHEM_ and QOCT-GRAD_. We also use it for the job scripts, containing all instructions for the job scheduler (SLURM_ for the CECI_ clusters).

Here is a small YouTube tutorial made by Jason Rigden that will cover the basics of the Jinja2 templating engine:

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/bxhXQG1qJPM" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

.. .. youtube:: bxhXQG1qJPM

|
| You can also find a written version of that video `here <https://medium.com/@jasonrigden/jinja2-templating-engine-tutorial-4bd31fb4aea3>`_. 

With that, you know everything you need to know about Jinja2 in order to use CHAINS! As you may have seen in the video, the Jinja language offers a relatively easy and intuitive way of creating those templates, which makes it easier to adapt the templates to pretty much any input files for any ab initio program such as ORCA_ or Molpro_, and any other job scheduler such as Torque_.

.. _regex:

About Regular Expressions
=========================

`Regular Expressions`_ (or regexes) are a powerful tool that can define search patterns to look for specific information in text. You don't need to know how they work to simply use ``ABIN LAUNCHER`` or ``CONTROL LAUNCHER``, but they are used in the results treatment part of CHAINS. If you want to customize the treatment of your results, or if you want to define new :doc:`scanning functions <abin_launcher.scan>` for geometry files (in ``ABIN LAUNCHER``) or new :ref:`modelling functions <modelling_fcts>` for source files (in ``CONTROL LAUNCHER``), it is probably a good idea to take a look at how they work.

Regular Expressions can present themselves as incredibly unintuitive but their versatility certainly makes up for that. Here is a YouTube tutorial made by Corey Schafer that will cover how regexes work and how they can be used in Python:

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/K8L6KVGG-7o" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

|
| With that, you know everything you need to know about regular expressions in order to customize CHAINS! To make it easier to work with regexes, feel free to use the `regex101 website <https://regex101.com/>`_ which can help you to build and to understand them.

.. _cron_tuto:

About Crontab and Cron tasks
============================

In order to link the different clusters and to automate certain parts of CHAINS, we need to periodically execute some commands. This is done through cron_ tasks established via the ``crontab`` command. 

You don't really need to know anything about them prior to reading this documentation, as the exact commands will be explicitly given each time you will need to use them. Nonetheless, if you prefer to know what you are dealing with, here is a small YouTube tutorial made by Corey Schafer explaining what cron tasks are:

.. raw:: html

   <iframe width="560" height="315" src="https://www.youtube.com/embed/QZJ1drMQz1A" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

|

.. Hyperlink targets

.. _CECI: http://www.ceci-hpc.be/
.. _Cron: https://pubs.opengroup.org/onlinepubs/9699919799/utilities/crontab.html
.. _Gaussian: https://gaussian.com/
.. _Jinja: https://jinja.palletsprojects.com/en/2.11.x/ 
.. _Molpro: https://www.molpro.net/
.. _ORCA: https://www.faccts.de/orca/
.. _Python: https://www.python.org/
.. _PyYAML: https://pyyaml.org/
.. _Q-CHEM: https://www.q-chem.com/
.. _QOCT-GRAD: https://gitlab.com/dynaq.cqp/QOCT-GRAD
.. _`Regular Expressions`: https://www.regular-expressions.info/
.. _SLURM: https://slurm.schedmd.com/documentation.html
.. _Torque: https://github.com/adaptivecomputing/torque
.. _YAML: https://yaml.org/