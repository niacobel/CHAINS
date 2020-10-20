.. CHAINS documentation master file, created by
   sphinx-quickstart on Sun Sep 27 16:57:36 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CHAINS' documentation!
=================================

CHAINS is a group of scripts developed at the Spectroscopy, Quantum Chemistry and Atmospheric Remote Sensing (SQUARES_) service, at *Université Libre de Bruxelles* (ULB_). Its use is to automate the scanning of multiple molecules, looking for the ones whose electrons can be easily controlled by laser. 

This includes two big steps: **caracterisation** and **evaluation**. The caracterisation step makes use of ab initio calculations performed by programs such as ORCA_ and Q-CHEM_, in order to build an effective Hamiltonian describing our molecule. The evaluation step consists of applying the quantum optimal control theory through the QOCT-RA_ scripts, also developed by SQUARES. All those calculations are performed using the different clusters provided by the *Consortium des Équipements de Calcul Intensif* (CECI_). 

Although the global development has been aimed at a specific problematic, **each individual script can easily be adapted** to deal with other similar problems, programs and clusters, as it will be explained in each of the different sections of this documentation.

Dependencies
------------

COMING SOON

Design Philosophy
-----------------

COMING SOON

Global Overview
---------------

.. figure:: figures/chains_workflow.*
    :width: 940px
    :align: center
    :height: 696px
    :alt: CHAINS workflow
    :target: _images/chains_workflow.png
    :figclass: align-center

    
    Global overview of CHAINS' workflow (click to enlarge)

Caracterisation Step (ABIN LAUNCHER)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step consists to scan the molecule structure files and building the input files associated with the ab initio program you want to run. Once those files have been prepared, the corresponding calculations will be launched. This whole process is covered by the first main script of CHAINS, called ABIN LAUNCHER. 

This script is executed twice through CHAINS. The first time, it is used to perform the geometry optimization of our molecule through the ORCA_ program. The second time, it is used to calculate the different properties we will need during the evaluation step, through the Q-CHEM_ program.

*Note: The reason this caracterisation step is split in two phases is because the Q-CHEM_ program is locked by license on a single cluster, and we want to make the most use out of this cluster. Since the geometry optimization can be handled by another program, we have separated it from the rest in order to free some ressources and gain some time.*

Evaluation Step (CONTROL LAUNCHER)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second step consists to scan the output file of the caracterisation step and use it to build an effective Hamiltonian for our molecule. Once in possession of this effective Hamiltonian, we then perform a quantum control procedure through QOCT-RA_, consisting of populating each dark electronic state, from the ground state. This whole process is covered by the second main script of CHAINS, called CONTROL LAUNCHER. Each of those calculations yields a fidelity value, which reflects the efficiency of the control procedure and the controllability of the electrons in our molecule.

Final Step (RESULTS TREATMENT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is no calculation involved in this step. The third main script of CHAINS, called RESULTS TREATMENT, is designed to, as its name implies, treat the results from the previous steps and compile them into graphs and tables, for ease of interpretation and comparison.

Linking the steps together
~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to link all those steps together and allow communication between clusters, CHAINS makes use of the common CECI storage, known as CECIHOME. Every important file is copied and stored into the CECIHOME, then different Shell scripts are periodically executed through crontab_ tasks to scan the CECIHOME and launch the various steps.

Getting Started
---------------

COMING SOON

Acknowledgment
--------------

COMING SOON

Copyright and License
---------------------

COMING SOON

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents

   Overview <cnt.overview>
   Configuration <cnt.config>
   Design Philosophy <cnt.design_phil>


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Docsite testing

   Tests <tests.rst>

.. _SQUARES: https://www2.ulb.ac.be/cpm/index.html
.. _ULB: https://www.ulb.be/
.. _ORCA: https://www.faccts.de/orca/
.. _Q-CHEM: https://www.q-chem.com/
.. _QOCT-RA: https://gitlab.com/dynaq.cqp/QOCT-RA
.. _CECI: http://www.ceci-hpc.be/
.. _crontab: https://pubs.opengroup.org/onlinepubs/9699919799/utilities/crontab.html