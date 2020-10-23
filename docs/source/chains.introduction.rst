************
Introduction
************

What is CHAINS?
===============

CHAINS is a group of scripts developed at the Spectroscopy, Quantum Chemistry and Atmospheric Remote Sensing (SQUARES_) service, at *Université Libre de Bruxelles* (ULB_). Its purpose is to automate a "quantum control screening" methodology, which consists of scanning multiple molecules, looking for the ones whose electrons can be easily controlled by laser.

This methodology includes three big phases:

- The **characterization** phase makes use of ab initio calculations performed by programs such as ORCA_ and Q-CHEM_, in order to build an effective Hamiltonian describing our molecule. 
- The **evaluation** phase consists of using this effective Hamiltonian to apply the quantum optimal control theory through the QOCT-RA_ scripts, also developed at SQUARES.
- The **results** phase consists of treating the results of all calculations and compiling them into graphs and tables.

CHAINS can assist this whole methodology by automatically preparing the input files and launching the ORCA_, Q-CHEM_ and QOCT-RA_ calculations for each molecule, allowing the user to treat a big number of molecules much faster, thus making the scanning more efficient. It can also work with multiple clusters, in order to split the workload. In our case, all those calculations (jobs) are performed on the different clusters provided by the *Consortium des Équipements de Calcul Intensif* (CECI_). 

.. note::
   Although the global development has been aimed at a specific problematic, **each individual script can easily be adapted** to deal with other similar problems, programs and clusters, as will be explained in each of the different sections of this documentation.

Design Philosophy
=================

The development of CHAINS followed two main ideas:

- **Simple**: Since CHAINS might be used by scientists with minimal knowledge of software development, we want to keep it simple and easy to understand, so we tried to stick with the basics as much as possible. Also, since the jobs that will be launched might last for hours, days and maybe even weeks, there is no need to put too much emphasis on efficiency, as gaining a few seconds during the execution of CHAINS' scripts wouldn't really impact the overall performance of the total procedure. 
- **Adaptable**: Since there are as many scientific problems as there are scientists, we want to make CHAINS easily adaptable. Therefore, we also tried to keep in mind other programs and methodologies that may want to use some or all of CHAINS' scripts. 

As a final note, CHAINS is not a single program but really just a group of scripts. Those scripts will be executed before and after each calculation, in order to prepare input files and treat output files. As such, if one calculation happens to fail, it is easy to restart where it stopped.

Languages and Dependencies
==========================

.. figure:: figures/logos.*
    :align: center
    :alt: Languages used in CHAINS
    :figclass: align-center

    The different languages used by CHAINS

The main language of CHAINS is Python 3.5+ but it also uses YAML files and Jinja2 templates, as well as some Shell scripts. However, there is no need for an in-depth knowledge of anything but Python in order to customize CHAINS, as we only use the basic features of the other languages.

In order to plot graphs and generate tables with the results, CHAINS also makes use of Gnuplot_ scripts and the LaTeX_ language. Those are not needed to get the results though, only to process them and make them more readable. You can freely rework those scripts and files without affecting the main process.

*What is YAML and why is it used here?*
---------------------------------------

YAML_ is a "human friendly data serialization standard for all programming languages". It is used here for the configuration files, containing all parameters needed for the different scripts and calculations. The way those files are structured makes them very easy to use and to understand.

*What is Jinja and why is it used here?*
----------------------------------------

"Jinja_ is a modern and designer-friendly templating language for Python, modelled after Django’s templates. It is fast, widely used and secure with the optional sandboxed template execution environment." (quote from the official website)

We use the Jinja language to create templates for the input files needed by ORCA_, Q-CHEM_ and QOCT-RA_. We also use it for the job scripts, containing all instructions for the job scheduler (SLURM_ for the CECI clusters). The Jinja language offers a relatively easy and intuitive way of creating those templates, which makes it easier to adapt the templates to pretty much any input files for any ab initio program such as Gaussian_ or Molpro_, and any other job scheduler such as Torque_.

*Installing YAML and Jinja2*
----------------------------

If your Python distribution does not include YAML and Jinja2, you can install them using

.. code-block:: shell

   $ python -m pip install pyyaml
   $ python -m pip install jinja2

Note that depending on your machine configuration, you might need to specify to install them only for your user account:

.. code-block:: shell

   $ python -m pip install --user pyyaml
   $ python -m pip install --user jinja2

Otherwise, you might be denied the permission to install.

Acknowledgment
==============

CHAINS makes use of a YAML version of Mendeleev's periodic table procured by `AlexGustafsson's molecular-data Github repository`_.

The main developer of CHAINS, Nicolas Iacobellis, would also like to express his deepest gratitude and give a shout-out to his friend, `Benjamin D'Heure`_, for its tremendous help and essential expertise during the code development. Without him, this project might not have existed, or would have at least taken a different form.

License
=======

.. todo::
   COMING SOON (Probably just GPLv3)

.. Hyperlink targets

.. _`AlexGustafsson's molecular-data Github repository`: https://github.com/AlexGustafsson/molecular-data
.. _`Benjamin D'Heure`: https://www.linkedin.com/in/bdheure/
.. _CECI: http://www.ceci-hpc.be/
.. _Gaussian: https://gaussian.com/
.. _Gnuplot: http://www.gnuplot.info/
.. _Jinja: https://jinja.palletsprojects.com/en/2.11.x/ 
.. _LaTeX: https://www.latex-project.org/
.. _Molpro: https://www.molpro.net/
.. _ORCA: https://www.faccts.de/orca/
.. _Q-CHEM: https://www.q-chem.com/
.. _QOCT-RA: https://gitlab.com/dynaq.cqp/QOCT-RA
.. _SLURM: https://slurm.schedmd.com/documentation.html
.. _SQUARES: https://www2.ulb.ac.be/cpm/index.html
.. _Torque: https://github.com/adaptivecomputing/torque
.. _ULB: https://www.ulb.be/
.. _YAML: https://yaml.org/