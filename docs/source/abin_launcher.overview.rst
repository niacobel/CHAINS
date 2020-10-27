*************************
Overview of ABIN LAUNCHER
*************************

What is ABIN LAUNCHER?
======================

``ABIN LAUNCHER`` is a script that creates input files for a given ab initio program, and then submits (or "launches") the corresponding calculation to a job scheduler. It can operate with one or more molecules (through their respective geometry files) and one or more configurations (the parameters for the ab initio program). ``ABIN LAUNCHER`` is the first main script of CHAINS and is executed twice through it:

- The first time, it is used to build the input files associated with the ORCA_ program, which performs the geometry optimization of the molecule. 
- The second time, it is used to build the input files associated with the the Q-CHEM_ program, which performs the calculation of the different properties of our molecule. 

In both cases, after having created the files, ``ABIN LAUNCHER`` then launches the corresponding job on the cluster.

Note however that ``ABIN LAUNCHER`` is **completely autonomous** and does not depend on CHAINS. It does not need any files outside the ones present in its own directory. It can be extracted and used to launch independent calculations, and can be very easily adapted to deal with other ab initio programs. For more information about how to adapt ``ABIN LAUNCHER`` to outside cases, consult the :doc:`Usage outside of CHAINS <abin_launcher.adapt>` specific documentation.

How does it work?
=================

``ABIN LAUNCHER``'s procedure follows three main steps: **scanning**, **scaling** and **rendering**.

Scanning
--------

``ABIN LAUNCHER`` begins by scanning the geometry file, looking for the chemical formula and the atomic coordinates of the molecule. For more details on how this scan is performed, consult the :doc:`Scanning the geometry file <abin_launcher.scan>` specific documentation.

.. note::
   At this time, only the XYZ format is supported for geometry files. However, new formats can be added in the future if the need arises.

Scaling
-------

Based on the information received from the geometry file, ``ABIN LAUNCHER`` attributes a value, called the scale index, to the molecule. This value is then used to evaluate the job scale for that molecule and specify the calculation requirements accordingly (walltime, number of CPUs, memory, etc.). For more details on how this scaling process is performed, consult the :doc:`Job scaling <abin_launcher.job_scale>` specific documentation.

Rendering
---------

Finally, based on user-defined Jinja templates, ``ABIN LAUNCHER`` creates the input files and the job instructions file associated with our calculation, then submits the calculation to the job scheduler. The content of those input files is based on the information from the geometry file and the configuration file (YAML file containing all the parameters for the ab initio calculation). For more details on how this whole rendering process is performed, consult the :doc:`Rendering the templates <abin_launcher.rendering>` specific documentation.

Directory structure
===================

The ``ABIN LAUNCHER`` directory has and must have the following structure:

.. code-block::

    abin_launcher/
      ├── abin_launcher.py
      ├── mol_scan.py
      ├── scaling_fcts.py
      ├── renderer.py
      ├── abin_errors.py
      ├── clusters.yml
      ├── mendeleev.yml
      └── Templates/
          ├── orca_job.sh.jinja
          ├── orca.inp.jinja
          ├── qchem_job.sh.jinja
          └── qchem.in.jinja

As for what each file does, everything will be explained in more details in the other sections of this documentation, but here is a short summary:

- ``abin_launcher.py`` is the main script itself, the one that needs to be executed.
- ``mol_scan.py`` is the library of functions that define how to scan the geometry files, i.e. how to read and interpret them.
- ``scaling_fcts`` is the library of functions that define how to determine the job size, by calculating what we defined as the scale index.
- ``renderer.py`` is the library of functions that define how to render the Jinja templates, i.e. to create the input files and the job instructions file.
- ``abin_errors.py`` contains all the classes and functions defining how to handle errors.
- ``clusters.yml`` is the YAML file containing all the information specific to the different clusters, notably the job scales and their associated requirements.
- ``mendeleev.yml`` is a YAML version of Mendeleev's periodic table procured by `AlexGustafsson's molecular-data Github repository`_

- Finally, ``Templates`` is the directory containing all the Jinja templates that will be used by ``ABIN LAUNCHER``

    - ``orca.inp.jinja`` and ``qchem.in.jinja`` are the Jinja templates corresponding to the input file for ORCA and Q-CHEM, respectively.
    - ``orca_job.sh.jinja`` and ``qchem_job.sh.jinja`` are the Jinja templates corresponding to the job instructions file for ORCA and Q-CHEM, respectively. They contain all the commands that will be executed through the job scheduler to perform the calculation.

.. warning::

   While you might be tempted to modify some of those scripts to adapt CHAINS to your case study, you should not alter ``mendeleev.yml`` as this serves as reference for key sections of the code.

Command line arguments
======================

.. argparse::
   :module: abin_launcher
   :func: parser
   :prog: abin_launcher.py
   :nodescription:

Usage
=====

.. todo::
   COMING SOON

.. Hyperlink targets

.. _`AlexGustafsson's molecular-data Github repository`: https://github.com/AlexGustafsson/molecular-data
.. _ORCA: https://www.faccts.de/orca/
.. _Q-CHEM: https://www.q-chem.com/
