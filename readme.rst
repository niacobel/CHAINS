.. badges

|GitHub Page| |Documentation Status| |GitHub License| |GitHub issues|

.. |GitHub Page| image:: https://img.shields.io/website-up-down-green-red/https/github.com/niacobel/CHAINS.svg
   :alt: GitHub Page
   :target: https://github.com/niacobel/CHAINS

.. |Documentation Status| image:: https://readthedocs.org/projects/chains-ulb/badge/
    :alt: Documentation Status
    :target: https://chains-ulb.readthedocs.io/en/latest/?badge=latest

.. |GitHub License| image:: https://img.shields.io/github/license/niacobel/CHAINS.svg
   :alt: GitHub License
   :target: https://github.com/niacobel/CHAINS/blob/master/LICENSE

.. |GitHub issues| image:: https://img.shields.io/github/issues/niacobel/CHAINS.svg
   :target: https://github.com/niacobel/CHAINS/issues/

What is CHAINS?
===============

CHAINS is a group of scripts developed at the Spectroscopy, Quantum Chemistry and Atmospheric Remote Sensing (SQUARES_) service, at *Université Libre de Bruxelles* (ULB_). Its purpose is to automate a "quantum control screening" methodology, which consists in scanning multiple molecules, looking for the ones whose electrons can be easily controlled by laser.

This methodology includes three big phases:

#. The **characterization** phase makes use of ab initio calculations performed by programs such as GAUSSIAN_ and Q-CHEM_, in order to build a matrix Hamiltonian describing our molecule. 
#. The **control** phase consists in using this Hamiltonian to apply the quantum optimal control theory through the QOCT-GRAD_ scripts, also developed at SQUARES.
#. The **results** phase involves treating the results of all calculations and compiling them into graphs and tables.

CHAINS can assist this whole methodology by automatically preparing the input files and launching the GAUSSIAN_, Q-CHEM_ and QOCT-GRAD_ calculations for each molecule, allowing the user to treat a big number of molecules much faster, thus making the scanning more efficient. It can also work with multiple clusters, in order to split the workload. In our case, all those calculations (jobs) are performed on the different clusters provided by the *Consortium des Équipements de Calcul Intensif* (CECI_). 

.. Important::
   Although the global development has been aimed at a specific problematic, **each individual script can easily be adapted** to deal with other similar problems, programs and clusters, as will be explained in each of the different sections of the documentation.

Design Philosophy
=================

The development of CHAINS followed two main ideas:

* **Simple**: Since CHAINS might be used by scientists with minimal knowledge of software development, we want to keep it simple and easy to understand, so we tried to stick with the basics as much as possible. Also, since the jobs that will be launched might last for hours, days and maybe even weeks, there is no need to put too much emphasis on efficiency, as gaining a few seconds during the execution of CHAINS' scripts wouldn't really impact the overall performance of the total procedure. 
* **Adaptable**: Since there are as many scientific problems as there are scientists, we want to make CHAINS easily adaptable. Therefore, we also tried to keep in mind other programs and methodologies that may want to use some or all of CHAINS' scripts. 

As a final note, CHAINS is not a single program but really just a group of scripts. Those scripts will be executed on the clusters before and after each calculation, in order to prepare input files and treat output files. As such, if one calculation happens to fail, it is easy to restart juste before that step.

Languages and Dependencies
==========================

.. figure:: https://raw.githubusercontent.com/niacobel/CHAINS/master/docs/source/figures/logos.png
    :align: center
    :alt: Languages used in CHAINS
    :figclass: align-center
    
|
The main language of CHAINS is **Python 3.5+** but it also uses YAML_ files and Jinja2_ templates, as well as some Bash scripts. If you want to customize CHAINS however, there is no need for an in-depth knowledge of anything but Python, as we only use the basic features of the other languages (consult the "What you need to know" section of the documentation for details).

If your Python distribution does not include them, you can install the needed Python libraries by entering the following command in your terminal:

.. code-block:: console

   $ python -m pip install --upgrade --user pyyaml jinja2 scipy

The third main script of CHAINS, ``RESULTS TREATMENT``, makes use of the Pandas_ and Matplotlib_ packages. If your Python distribution does not include them, you can install them using

.. code-block:: console

   $ python -m pip install --user pandas matplotlib

.. Important::
   Since CHAINS' scripts will be running on the clusters themselves and not on your personal computer, you need to make sure to meet those packages requirements for each cluster on which you intend to use CHAINS.

Acknowledgment
==============

CHAINS makes use of a YAML version of Mendeleev's periodic table procured by `AlexGustafsson's molecular-data Github repository`_.

The main developer of CHAINS, Nicolas Iacobellis, would also like to express his deepest gratitude and give a shout-out to his friend, Benjamin D'Heure (`Linkedin page`_), for its tremendous help and essential expertise during the code development. Without him, this project might not have existed, or would have at least taken a different form.

License
=======

.. todo::
   COMING SOON (Probably just GPLv3)

.. Hyperlink targets

.. _`AlexGustafsson's molecular-data Github repository`: https://github.com/AlexGustafsson/molecular-data
.. _`Linkedin page`: https://www.linkedin.com/in/bdheure/
.. _CECI: http://www.ceci-hpc.be/
.. _GAUSSIAN: https://gaussian.com/
.. _Jinja2: https://jinja.palletsprojects.com/en/2.11.x/ 
.. _LaTeX: https://www.latex-project.org/
.. _Matplotlib: https://matplotlib.org/
.. _Pandas: https://pandas.pydata.org/
.. _Q-CHEM: https://www.q-chem.com/
.. _QOCT-GRAD: https://gitlab.com/dynaq.cqp/QOCT-GRAD
.. _SciPy: https://scipy.org/
.. _SQUARES: https://www2.ulb.ac.be/cpm/index.html
.. _ULB: https://www.ulb.be/
.. _YAML: https://yaml.org/