********
Workflow
********

.. figure:: figures/workflow.*
    :scale: 65%
    :align: center
    :alt: CHAINS workflow
    :figclass: align-center

    
    Global overview of CHAINS' workflow (click to zoom in)

Characterization Phase (Steps 1-3 and 5-6)
==========================================

**This phase is covered by the first main script of CHAINS, named** ``ABIN LAUNCHER``.

``ABIN LAUNCHER`` starts by scanning the molecule structure file (the geometry file) and building the input files associated with the ab initio program we want to run. Once those files have been prepared, the corresponding job will be launched. Note that the characterization phase is split into two jobs and that ``ABIN LAUNCHER`` is executed twice. The first time, it is used to perform the geometry optimization of the molecule through the ORCA_ program, and the second time, it is used to calculate different properties of our molecule through the Q-CHEM_ program.

.. note::
   ``ABIN LAUNCHER`` **is completely autonomous.** It does not need any files outside the ones present in its own directory. It can be extracted and used to launch independent calculations, and can be very easily adapted to deal with other ab initio programs. Consult ``ABIN LAUNCHER``'s specific documentation for more information.

.. note:: 
   The reason this characterization phase is split in two jobs is because the Q-CHEM_ program is locked by license on a single cluster, and we want to make the most use out of this cluster. Since the geometry optimization can be handled by another program, we have separated it from the rest in order to free some resources and gain some time.

Evaluation Phase (Steps 8-9)
============================

**This phase is covered by the second main script of CHAINS, named** ``CONTROL LAUNCHER``.

``CONTROL LAUNCHER`` starts by scanning the Q-CHEM output file and use it to build an effective Hamiltonian for our molecule. Once in possession of this effective Hamiltonian, multiple jobs will be launched, in order to perform a quantum control procedure through QOCT-RA_, consisting of separately populating each dark electronic state, starting from the ground state. Each of those calculations yields a fidelity value, which reflects the efficiency of the control procedure and the controllability of the electrons in our molecule.

Results Phase (Step 11)
=======================

**This phase is covered by the third main script of CHAINS, named** ``RESULTS TREATMENT``.

There is no calculation involved in this phase. As its name implies, ``RESULTS TREATMENT`` is designed to treat the results from the previous steps and compile them into graphs and tables, for ease of interpretation and comparison.

Link between scripts (Steps 4, 7 and 10)
========================================

In order to link all those scripts together and allow communication between clusters, CHAINS makes use of the common CECI storage, known as ``CECIHOME``. Every important file is copied and stored into the ``CECIHOME``, then different Shell scripts on different clusters are periodically executed through crontab_ tasks to scan the ``CECIHOME`` and launch the various tasks.

.. Hyperlink targets

.. _crontab: https://pubs.opengroup.org/onlinepubs/9699919799/utilities/crontab.html
.. _ORCA: https://www.faccts.de/orca/
.. _Q-CHEM: https://www.q-chem.com/
.. _QOCT-RA: https://gitlab.com/dynaq.cqp/QOCT-RA
