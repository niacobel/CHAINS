****************
System modelling
****************

``CONTROL LAUNCHER`` begins by modelling the system on which the control procedure will be performed. In order to do so, it must parse the source file to extract all the needed values from it and build the Matrix Image of the MoleculE (MIME), as well as the transition dipole moments matrix. The MIME can act as an effective Hamiltonian describing the molecule. Once it has been determined, it is then diagonalized in order to build the eigenstates basis set needed by QOCT-RA to perform the control procedure.

.. _parsing_fcts:

Parsing the source file
=======================

While ``CONTROL LAUNCHER`` parses the source file, it is looking for various information about the molecule, depending on what type of control is to be performed. What this parsing does and how it works is defined by a function present in the ``source_parser.py`` file, called a **parsing function**. The only argument this function needs is the content of the source file, stored as a list in the ``source_content`` variable (each element of the list corresponds to a line of the file). 

General definition of the parsing functions
-------------------------------------------

All the parsing functions are defined in the ``source_parser.py`` file and obey some restrictions, in order to be callable by ``CONTROL LAUNCHER``:

- They must only take one argument: a list containing the lines of the source file (``source_content``).
- They must return a dictionary (``system``) containing two mandatory keys: ``mime`` and ``momdip_mtx`` (you can have additional keys if you want) where:

   - ``mime`` is the Matrix Image of the MoleculE, acting as an effective Hamiltonian. The diagonal elements are the energy values of the states and the non-diagonal values are the coupling elements between those states. Those values must be expressed in **atomic units of energy (Hartree)**.
   - ``momdip_mtx`` is the transition dipole moments matrix containing the transition dipole moments between the states. Those values must be expressed in **atomic units**.

The values for ``system``'s two keys can either be NumPy 2D-arrays or a simple list of lists, depending on your preferences.

If a problem arises when parsing the source file, a ``ControlError`` exception should be raised with a proper error message (see :ref:`how to handle errors <control_errors>` for more details).

Choosing a parsing function
---------------------------

The parsing function that will be called by ``CONTROL LAUNCHER`` is the one associated with the ``parsing_function`` YAML key defined in the :ref:`clusters configuration file <control_clusters_file>`:

.. code-block:: yaml

   mycluster:
      profiles:
         myprofile1:
            parsing_function: name-of-parsing-function
         myprofile2:
            parsing_function: name-of-parsing-function

where ``mycluster`` corresponds to the name of your cluster (given as a :ref:`command line argument <control_arguments>`) while ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to run (such as ``chains_qoctra``). This way, a different parsing function can be assigned to each profile.

Convert energy units
--------------------

Since the ``mime`` values need to be expressed in Hartree, a small function named ``energy_unit_conversion`` has been added at the beginning of the ``source_parser.py`` file to help you convert energy units:

.. autofunction:: source_parser.energy_unit_conversion

This function works by first converting the value from the ``init`` unit to Hartree, then converting it from Hartree to the ``target`` unit.

As an example, let's say we want to convert 50 eV into Hz, we can simply use:

.. code-block:: python

   value_hz = energy_unit_conversion(50,"eV","Hz")

which gives us ~1.209e+16 Hz. If you simply want the value in Hartree, just use ``"Ha"`` as the target unit (which would give ~1.837 Ha in this case).

Note that the unit labels are case insensitive. If you want to add support for another unit, simply add the conversion factor from Hartree to that unit to the ``conv_factors`` dictionary defined in the function.

Q-CHEM TDDFT parsing function
-----------------------------

This parsing function scans the content of a Q-CHEM TD-DFT calculation output file. Note that the ground state of the calculation must be a singlet, otherwise this function will need to be slightly modified.

Structure of the Q-CHEM output file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, let's consider here a TD-DFT calculation of the |Si5OH10| molecule, focused on four excited triplet states and four excited singlet states. The output file is quite lengthy so we will only focus on the relevant parts in this subsection but know that the whole file can be consulted `here <https://github.com/niacobel/CHAINS/tree/master/docs/source/control_sample/sp-si5oh10.out>`_.

.. |Si5OH10| replace:: Si\ :sub:`5`\ OH\ :sub:`10`\ 

.. container:: toggle

   .. container:: header

      The first part of interest is the list of the **excited electronic states and their excitation energies**. *(click on the arrow below to see this portion of the output file.)*

   .. literalinclude:: control_sample/si5oh10.out
      :language: text
      :caption: si5oh10.out - TDDFT/TDA Excitation Energies
      :lines: 276-344

|

.. container:: toggle

   .. container:: header

      Next comes the **spin-orbit couplings** between those excited states (and the ground state). *(click on the arrow below to see this portion of the output file.)*

   .. literalinclude:: control_sample/si5oh10.out
      :language: text
      :caption: si5oh10.out - Spin-orbit Couplings
      :lines: 348-558

|

.. container:: toggle

   .. container:: header

      And finally, the **transition dipole moments** between those electronic states. *(click on the arrow below to see this portion of the output file.)*

   .. literalinclude:: control_sample/si5oh10.out
      :language: text
      :caption: si5oh10.out - Transition dipole moments
      :lines: 563-643

|

Function definition
~~~~~~~~~~~~~~~~~~~

.. autofunction:: source_parser.qchem_tddft

The first noticeable thing is that the ``system`` dictionary returned by this function contains an additional key: ``states_list``. While this key is not explicitly needed by ``CONTROL LAUNCHER``, it contains useful information that can be used, for example, by a transition function.

The way this function works is by making heavy use of :ref:`regular expressions <regex>`. It detects the beginning and the end of each of the three portions presented above and process their content, picking the relevant values and converting them to atomic units if needed.

Example with the source file shown in the previous subsection:

.. code-block:: text

   Parsing the excited states ...                   [ DONE ]

   --------------------------------------------------
                     States List
   --------------------------------------------------
   Number     Multiplicity    Energy (cm-1)   Label
   --------------------------------------------------
   0          Singlet         0.000           S0
   1          Triplet         37042.624       T1
   2          Triplet         38400.861       T2
   3          Triplet         38692.027       T3
   4          Triplet         42177.149       T4
   5          Singlet         42192.473       S1
   6          Singlet         42956.280       S2
   7          Singlet         45173.498       S3
   8          Singlet         46441.402       S4
   --------------------------------------------------

   Parsing the spin-orbit couplings ...             [ DONE ]

   Building the MIME ...                            [ DONE ]

   MIME (cm-1)

   0.00000e+00  4.18679e+01  1.60967e+00  1.70877e+01  1.98183e+01  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
   4.18679e+01  3.70426e+04  6.04852e+01  2.64589e+01  2.40051e+01  1.87819e-01  1.72586e+01  5.64945e+01  2.37051e+01
   1.60967e+00  6.04852e+01  3.84009e+04  3.91753e+01  9.74640e+00  1.30716e+01  2.06946e+01  1.01841e+00  2.66124e+00
   1.70877e+01  2.64589e+01  3.91753e+01  3.86920e+04  7.26096e+01  1.28708e+01  9.43968e-01  2.54447e+00  1.33685e+01
   1.98183e+01  2.40051e+01  9.74640e+00  7.26096e+01  4.21771e+04  1.78498e+01  5.36720e+01  1.89108e+01  2.19255e+01
   0.00000e+00  1.87819e-01  1.30716e+01  1.28708e+01  1.78498e+01  4.21925e+04  0.00000e+00  0.00000e+00  0.00000e+00
   0.00000e+00  1.72586e+01  2.06946e+01  9.43968e-01  5.36720e+01  0.00000e+00  4.29563e+04  0.00000e+00  0.00000e+00
   0.00000e+00  5.64945e+01  1.01841e+00  2.54447e+00  1.89108e+01  0.00000e+00  0.00000e+00  4.51735e+04  0.00000e+00
   0.00000e+00  2.37051e+01  2.66124e+00  1.33685e+01  2.19255e+01  0.00000e+00  0.00000e+00  0.00000e+00  4.64414e+04

   MIME (Ha)

   0.00000e+00  1.90764e-04  7.33421e-06  7.78573e-05  9.02986e-05  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
   1.90764e-04  1.68779e-01  2.75591e-04  1.20556e-04  1.09375e-04  8.55766e-07  7.86362e-05  2.57408e-04  1.08008e-04
   7.33421e-06  2.75591e-04  1.74967e-01  1.78496e-04  4.44078e-05  5.95587e-05  9.42913e-05  4.64022e-06  1.21255e-05
   7.78573e-05  1.20556e-04  1.78496e-04  1.76294e-01  3.30834e-04  5.86436e-05  4.30103e-06  1.15934e-05  6.09115e-05
   9.02986e-05  1.09375e-04  4.44078e-05  3.30834e-04  1.92173e-01  8.13295e-05  2.44548e-04  8.61640e-05  9.99000e-05
   0.00000e+00  8.55766e-07  5.95587e-05  5.86436e-05  8.13295e-05  1.92243e-01  0.00000e+00  0.00000e+00  0.00000e+00
   0.00000e+00  7.86362e-05  9.42913e-05  4.30103e-06  2.44548e-04  0.00000e+00  1.95723e-01  0.00000e+00  0.00000e+00
   0.00000e+00  2.57408e-04  4.64022e-06  1.15934e-05  8.61640e-05  0.00000e+00  0.00000e+00  2.05826e-01  0.00000e+00
   0.00000e+00  1.08008e-04  1.21255e-05  6.09115e-05  9.99000e-05  0.00000e+00  0.00000e+00  0.00000e+00  2.11603e-01

   Parsing the transition dipole moments ...        [ DONE ]

   Dipole moments matrix (atomic units)

   0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  5.56471e-02  5.73960e-05  6.91184e-04  6.53991e-02
   0.00000e+00  0.00000e+00  6.79321e-05  2.28178e-05  3.70329e-04  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
   0.00000e+00  6.79321e-05  0.00000e+00  4.74639e-07  3.66298e-04  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
   0.00000e+00  2.28178e-05  4.74639e-07  0.00000e+00  3.32524e-03  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
   0.00000e+00  3.70329e-04  3.66298e-04  3.32524e-03  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
   5.56471e-02  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  1.53868e-05  4.28151e-03  8.41199e-04
   5.73960e-05  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  1.53868e-05  0.00000e+00  2.67246e-07  1.38369e-05
   6.91184e-04  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  4.28151e-03  2.67246e-07  0.00000e+00  1.62057e-05
   6.53991e-02  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  8.41199e-04  1.38369e-05  1.62057e-05  0.00000e+00

Eigenstates basis set
=====================

.. todo::

   COMING SOON