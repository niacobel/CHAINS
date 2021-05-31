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

The values for ``system``'s two keys can either be NumPy_ 2D-arrays or a simple list of lists, depending on your preferences.

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

where ``mycluster`` corresponds to the name of your cluster (given as a :ref:`command line argument <control_arguments>`) while ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to use (such as ``chains_qoctra``). This way, a different parsing function can be assigned to each profile.

Convert energy units
--------------------

Since the ``mime`` values need to be expressed in Hartree, a small function named ``energy_unit_conversion`` has been added at the beginning of the ``source_parser.py`` file to help you convert energy units:

.. autofunction:: source_parser.energy_unit_conversion

This function works by first converting the value from the ``init`` unit to Hartree, then converting it from Hartree to the ``target`` unit.

As an example, let's say we want to convert 50 eV into Hz, we can simply use:

.. code-block:: python

   value_hz = energy_unit_conversion(50,"eV","Hz")

which gives us ~1.209e+16 Hz. If you simply want the value in Hartree, just use ``"Ha"`` as the target unit (which would give ~1.837 Ha in this case).

If you want to add support for another unit, simply add the conversion factor from Hartree to that unit to the ``conv_factors`` dictionary defined in the function. Note that the unit labels are case insensitive.

Q-CHEM TDDFT parsing function
-----------------------------

This parsing function uses a Q-CHEM_ TD-DFT calculation output as the source file. Note that the ground state of the calculation must be a singlet, otherwise this function will need to be slightly modified.

Structure of the Q-CHEM output file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, let's consider here a TD-DFT calculation of the |H2O| molecule, focused on four excited triplet states and four excited singlet states. The output file is quite lengthy so we will only focus on the relevant parts in this subsection, where the lines containing the needed information are emphasized. The whole file can however be consulted `here <https://github.com/niacobel/CHAINS/tree/master/docs/source/control_sample/h2o.out>`_.

.. |H2O| replace:: H\ :sub:`2`\ O

.. container:: toggle

   .. container:: header

      The first part of interest is the list of the **excited electronic states and their excitation energies**. *(click on the arrow below to see this portion of the output file.)*

   .. literalinclude:: control_sample/h2o.out
      :language: text
      :caption: h2o.out - TDDFT/TDA Excitation Energies
      :lines: 225-292
      :emphasize-lines: 5,7,12,14,19,21,26,28,33,35,40,42,48,50,56,58

|

.. container:: toggle

   .. container:: header

      Next comes the **spin-orbit couplings** between those excited states (and the ground state). *(click on the arrow below to see this portion of the output file.)*

   .. literalinclude:: control_sample/h2o.out
      :language: text
      :caption: h2o.out - Spin-orbit Couplings
      :lines: 296-506
      :emphasize-lines: 23-27,60-63,85-87,103-104,131-135,154-158,177-181,200-204

|

.. container:: toggle

   .. container:: header

      And finally, the **transition dipole moments** between those electronic states. *(click on the arrow below to see this portion of the output file.)*

   .. literalinclude:: control_sample/h2o.out
      :language: text
      :caption: h2o.out - Transition dipole moments
      :lines: 511-591
      :emphasize-lines: 27-30,37-42,61-64,71-76

|

Function definition
~~~~~~~~~~~~~~~~~~~

.. autofunction:: source_parser.qchem_tddft

The first noticeable thing is that the ``system`` dictionary returned by this function contains an additional key: ``states_list``. While this key is not explicitly needed by ``CONTROL LAUNCHER``, it contains useful information that can be used, for example, by a :ref:`transition function <pgt_transition>`.

The way this function works is by making heavy use of :ref:`regular expressions <regex>`. It detects the beginning and the end of each of the three portions presented above and process their content, picking the relevant values and converting them to atomic units if needed.

Example with the |H2O| output file:

.. code-block:: text

   Parsing the excited states ...                   [ DONE ]

   --------------------------------------------------
                     States List
   --------------------------------------------------
   Number     Multiplicity    Energy (cm-1)   Label
   --------------------------------------------------
   0          Singlet         0.000           S0
   1          Triplet         41673.859       T1
   2          Singlet         49010.278       S1
   3          Triplet         56123.281       T2
   4          Triplet         57236.326       T3
   5          Singlet         63149.176       S2
   6          Triplet         69737.113       T4
   7          Singlet         72968.976       S3
   8          Singlet         85210.859       S4
   --------------------------------------------------

   Parsing the spin-orbit couplings ...             [ DONE ]

   Building the MIME ...                            [ DONE ]

   MIME (cm-1)

   0.00000e+00  8.23059e+01  0.00000e+00  1.07552e+01  1.00613e+02  0.00000e+00  3.98429e+01  0.00000e+00  0.00000e+00
   8.23059e+01  4.16739e+04  2.42751e-01  1.12773e+02  6.14058e+01  4.21661e+01  5.13109e+01  7.88919e+01  4.10644e+01
   0.00000e+00  2.42751e-01  4.90103e+04  7.99241e+01  4.22015e+01  0.00000e+00  3.68396e+01  0.00000e+00  0.00000e+00
   1.07552e+01  1.12773e+02  7.99241e+01  5.61233e+04  1.47020e+01  9.83572e+00  1.76701e+01  1.54392e+01  1.96873e+00
   1.00613e+02  6.14058e+01  4.22015e+01  1.47020e+01  5.72363e+04  1.81588e-01  9.94394e+01  1.88645e+01  6.73752e+01
   0.00000e+00  4.21661e+01  0.00000e+00  9.83572e+00  1.81588e-01  6.31492e+04  7.07071e+01  0.00000e+00  0.00000e+00
   3.98429e+01  5.13109e+01  3.68396e+01  1.76701e+01  9.94394e+01  7.07071e+01  6.97371e+04  8.13273e+01  8.72452e+00
   0.00000e+00  7.88919e+01  0.00000e+00  1.54392e+01  1.88645e+01  0.00000e+00  8.13273e+01  7.29690e+04  0.00000e+00
   0.00000e+00  4.10644e+01  0.00000e+00  1.96873e+00  6.73752e+01  0.00000e+00  8.72452e+00  0.00000e+00  8.52109e+04

   MIME (Ha)

   0.00000e+00  3.75013e-04  0.00000e+00  4.90045e-05  4.58425e-04  0.00000e+00  1.81538e-04  0.00000e+00  0.00000e+00
   3.75013e-04  1.89880e-01  1.10605e-06  5.13834e-04  2.79785e-04  1.92123e-04  2.33790e-04  3.59458e-04  1.87103e-04
   0.00000e+00  1.10605e-06  2.23307e-01  3.64161e-04  1.92284e-04  0.00000e+00  1.67854e-04  0.00000e+00  0.00000e+00
   4.90045e-05  5.13834e-04  3.64161e-04  2.55716e-01  6.69873e-05  4.48148e-05  8.05110e-05  7.03464e-05  8.97019e-06
   4.58425e-04  2.79785e-04  1.92284e-04  6.69873e-05  2.60788e-01  8.27376e-07  4.53079e-04  8.59529e-05  3.06984e-04
   0.00000e+00  1.92123e-04  0.00000e+00  4.48148e-05  8.27376e-07  2.87729e-01  3.22165e-04  0.00000e+00  0.00000e+00
   1.81538e-04  2.33790e-04  1.67854e-04  8.05110e-05  4.53079e-04  3.22165e-04  3.17746e-01  3.70555e-04  3.97519e-05
   0.00000e+00  3.59458e-04  0.00000e+00  7.03464e-05  8.59529e-05  0.00000e+00  3.70555e-04  3.32471e-01  0.00000e+00
   0.00000e+00  1.87103e-04  0.00000e+00  8.97019e-06  3.06984e-04  0.00000e+00  3.97519e-05  0.00000e+00  3.88249e-01

   Parsing the transition dipole moments ...        [ DONE ]

   Dipole moments matrix (atomic units)

   0.00000e+00  0.00000e+00  3.76780e-03  0.00000e+00  0.00000e+00  2.86808e-04  0.00000e+00  1.01277e-01  1.33564e-02
   0.00000e+00  0.00000e+00  0.00000e+00  8.82015e-04  1.02219e-01  0.00000e+00  2.89492e-05  0.00000e+00  0.00000e+00
   3.76780e-03  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  9.64694e-02  0.00000e+00  1.96567e-03  9.02087e-07
   0.00000e+00  8.82015e-04  0.00000e+00  0.00000e+00  9.51555e-08  0.00000e+00  5.07242e-02  0.00000e+00  0.00000e+00
   0.00000e+00  1.02219e-01  0.00000e+00  9.51555e-08  0.00000e+00  0.00000e+00  7.04260e-04  0.00000e+00  0.00000e+00
   2.86808e-04  0.00000e+00  9.64694e-02  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00  6.53434e-07  8.57262e-04
   0.00000e+00  2.89492e-05  0.00000e+00  5.07242e-02  7.04260e-04  0.00000e+00  0.00000e+00  0.00000e+00  0.00000e+00
   1.01277e-01  0.00000e+00  1.96567e-03  0.00000e+00  0.00000e+00  6.53434e-07  0.00000e+00  0.00000e+00  1.07618e-01
   1.33564e-02  0.00000e+00  9.02087e-07  0.00000e+00  0.00000e+00  8.57262e-04  0.00000e+00  1.07618e-01  0.00000e+00

.. Caution::

   This function was used on output files generated by Q-CHEM versions 5.2.1 and 5.3. If you are using another version of Q-CHEM, the output might be different and the function might need to be adapted.

Eigenstates basis set
=====================

Once the MIME and the transition dipole moments matrix have been recovered from the source file, it is time to build the eigenstates basis set that will be used by QOCT-RA. In order to do so, the MIME needs to be diagonalized, this is done with the NumPy_ package, through the linalg.eig_ function.

.. code-block:: python

   system['eigenvalues'], system['eigenvectors'] = np.linalg.eig(system['mime'])

This gives us the eigenvalues, which are the eigenstates energies, and the associated eigenvectors, which are regrouped into a 2D-array. Those values are then sorted by ascending order of eigenvalues, so that the states are organized in ascending order of energies.

The last step is to compute the transpose of the eigenvectors matrix, with the transpose_ function from NumPy_:

.. code-block:: python

   system['transpose'] = np.transpose(system['eigenvectors'])

With the eigenvectors matrix and its transpose, we can now convert the transition dipole moments matrix to the eigenstates basis set, using the matmul_ function from NumPy_:

.. code-block:: python

   system['momdip_es_mtx'] = np.matmul(np.matmul(system['transpose'],system['momdip_mtx']),system['eigenvectors'])

Even though it's not needed, the diagonalized MIME is also recomputed, so that the diagonalization efficiency can be evaluated:

.. code-block:: python

   system['mime_diag'] = np.matmul(np.matmul(system['transpose'],system['mime']),system['eigenvectors'])

Note that all those computed variables have been added to the ``system`` dictionary, for ease of access by the other parts of ``CONTROL LAUNCHER``.

Let's end this subsection by using the |H2O| Q-CHEM_ output file as an example:

.. code-block:: text

   Diagonalizing the MIME ...                       [ DONE ]

   Eigenvalues (Ha)

   -1.6557807727135554e-06
   0.18987382916367557
   0.22330190065149785
   0.25572335800927826
   0.2607872862126559
   0.28772579137139187
   0.31774427201646704
   0.33248157890059926
   0.3882501806655301

   Eigenvectors matrix

   -9.99996e-01 -1.96175e-03 -1.46461e-05  1.77022e-04  1.76039e-03 -4.20907e-06  5.84862e-04 -1.87018e-05  4.04557e-06
    1.97120e-03 -9.99953e-01 -1.93928e-04  7.73046e-03  4.01626e-03  1.94498e-03  1.79630e-03 -2.56768e-03  9.47513e-04
   -2.25213e-06 -8.34760e-05  9.99923e-01  1.11428e-02  5.22847e-03 -2.04159e-05  1.79791e-03 -4.44360e-05  3.57291e-06
    1.87040e-04  7.79518e-03 -1.12157e-02  9.99808e-01  1.38601e-02  1.40381e-03  1.31095e-03 -9.61265e-04  7.29489e-05
    1.75469e-03  3.93215e-03 -5.08649e-03 -1.39604e-02  9.99845e-01 -1.26364e-04  7.93151e-03 -1.36897e-03  2.41264e-03
   -1.98561e-06  1.95372e-03  1.71601e-05 -1.43349e-03  1.32757e-05  9.99939e-01  1.07426e-02 -1.93645e-04  3.73071e-06
    5.67332e-04  1.79979e-03 -1.74281e-03 -1.24264e-03 -8.00107e-03 -1.07488e-02  9.99588e-01 -2.52312e-02  5.82656e-04
   -3.25671e-06  2.50978e-03  1.77862e-05 -9.30988e-04 -1.19125e-03  7.14264e-05 -2.52479e-02 -9.99677e-01  1.37866e-05
   -2.39976e-06  9.36334e-04  1.07164e-05 -4.58767e-05 -2.41244e-03  8.91066e-07 -6.03052e-04  3.42902e-05  9.99996e-01

   Eigenvectors transpose matrix

   -9.99996e-01  1.97120e-03 -2.25213e-06  1.87040e-04  1.75469e-03 -1.98561e-06  5.67332e-04 -3.25671e-06 -2.39976e-06
   -1.96175e-03 -9.99953e-01 -8.34760e-05  7.79518e-03  3.93215e-03  1.95372e-03  1.79979e-03  2.50978e-03  9.36334e-04
   -1.46461e-05 -1.93928e-04  9.99923e-01 -1.12157e-02 -5.08649e-03  1.71601e-05 -1.74281e-03  1.77862e-05  1.07164e-05
    1.77022e-04  7.73046e-03  1.11428e-02  9.99808e-01 -1.39604e-02 -1.43349e-03 -1.24264e-03 -9.30988e-04 -4.58767e-05
    1.76039e-03  4.01626e-03  5.22847e-03  1.38601e-02  9.99845e-01  1.32757e-05 -8.00107e-03 -1.19125e-03 -2.41244e-03
   -4.20907e-06  1.94498e-03 -2.04159e-05  1.40381e-03 -1.26364e-04  9.99939e-01 -1.07488e-02  7.14264e-05  8.91066e-07
    5.84862e-04  1.79630e-03  1.79791e-03  1.31095e-03  7.93151e-03  1.07426e-02  9.99588e-01 -2.52479e-02 -6.03052e-04
   -1.87018e-05 -2.56768e-03 -4.44360e-05 -9.61265e-04 -1.36897e-03 -1.93645e-04 -2.52312e-02 -9.99677e-01  3.42902e-05
    4.04557e-06  9.47513e-04  3.57291e-06  7.29489e-05  2.41264e-03  3.73071e-06  5.82656e-04  1.37866e-05  9.99996e-01

   MIME in the eigenstates basis set (Ha)

   -1.65578e-06  5.17176e-19  8.26105e-20 -3.83673e-20  4.84987e-20  1.62193e-22 -1.12580e-19  4.18511e-21 -9.63496e-23
    4.51967e-19  1.89874e-01 -1.88871e-16  2.14489e-18 -2.17604e-16 -1.00210e-19  2.86313e-17  4.52254e-18  6.06340e-17
    8.24916e-20 -1.88867e-16  2.23302e-01  1.07060e-16 -2.36417e-16  1.89052e-16 -9.51223e-17 -3.23605e-16  6.92267e-17
   -4.26371e-20  1.75558e-18  1.08101e-16  2.55723e-01  3.82438e-17 -5.42028e-16  4.64968e-17 -2.49204e-16  8.48270e-17
    1.37265e-19 -2.17695e-16 -2.36238e-16  3.71988e-17  2.60787e-01  6.43149e-17  7.40554e-16 -7.01983e-17  1.21322e-16
    4.63403e-22 -1.34276e-19  1.89053e-16 -5.42005e-16  6.43137e-17  2.87726e-01 -6.63569e-20 -1.31016e-16 -1.06336e-16
   -1.98754e-19  2.84508e-17 -9.50175e-17  4.66660e-17  7.40323e-16  3.34245e-20  3.17744e-01 -5.79015e-16 -2.00848e-17
    5.04310e-21  4.44232e-18 -3.23606e-16 -2.49150e-16 -7.00717e-17 -1.31017e-16 -5.77398e-16  3.32482e-01  6.61872e-18
   -2.11758e-22  6.06611e-17  6.92263e-17  8.48253e-17  1.21322e-16 -1.06336e-16 -2.00577e-17  6.62041e-18  3.88250e-01

   Dipole moments matrix in the eigenstates basis set (atomic units)

    1.46187e-06 -4.45419e-04 -3.77106e-03  8.23924e-05  3.36084e-04 -2.93927e-04  2.56796e-03  1.01243e-01 -1.33574e-02
   -4.45419e-04 -8.16837e-04  7.14051e-04  6.41318e-04 -1.02209e-01  8.10129e-07 -4.39198e-04  2.28362e-04 -4.88745e-07
   -3.77106e-03  7.14051e-04  5.47117e-06 -2.31070e-04 -1.42461e-05  9.64612e-02  4.15342e-04 -1.96738e-03  2.14454e-06
    8.23924e-05  6.41318e-04 -2.31070e-04 -1.37580e-04  3.80026e-04  5.28634e-04  5.06996e-02 -1.31406e-03 -6.81140e-05
    3.36084e-04 -1.02209e-01 -1.42461e-05  3.80026e-04  7.98671e-04  6.85852e-04  1.60022e-03 -2.27197e-04 -6.05293e-06
   -2.93927e-04  8.10129e-07  9.64612e-02  5.28634e-04  6.85852e-04 -5.51515e-06  2.44983e-04 -6.04909e-06  8.65636e-04
    2.56796e-03 -4.39198e-04  4.15342e-04  5.06996e-02  1.60022e-03  2.44983e-04  1.50952e-04 -5.19329e-05 -2.69342e-03
    1.01243e-01  2.28362e-04 -1.96738e-03 -1.31406e-03 -2.27197e-04 -6.04909e-06 -5.19329e-05 -1.78773e-07 -1.07585e-01
   -1.33574e-02 -4.88745e-07  2.14454e-06 -6.81140e-05 -6.05293e-06  8.65636e-04 -2.69342e-03 -1.07585e-01  3.55564e-06

Creating the files
==================

The last step of this system modelling process is to create the actual data files, writing down everything that has been computed until now. There are six of those files:

- ``mime`` contains the MIME in Hartree, stored in the ``mime`` key of the ``system`` dictionary.
- ``momdip_mtx`` contains the transition dipole moments matrix in atomic units, stored in the ``momdip_mtx`` key of the ``system`` dictionary.
- ``eigenvalues`` contains the list of the eigenvalues in Hartree, stored in the ``eigenvalues`` key of the ``system`` dictionary.
- ``eigenvectors`` contains the eigenvectors matrix, stored in the ``eigenvectors`` key of the ``system`` dictionary.
- ``transpose`` contains the eigenvectors transpose matrix, stored in the ``transpose`` key of the ``system`` dictionary.
- ``momdip_es_mtx`` contains the transition dipole moments matrix in atomic units, converted to the eigenstates basis set, stored in the ``momdip_es_mtx`` key of the ``system`` dictionary.

Those files are all created inside a new directory named ``data``, inside another new molecule directory bearing the name of the source file (minus a possible extension). At the end of the system modelling process, the output directory structure is then

.. code-block:: text

    out_dir/ 
      └── source/
            └── data/
                  ├── mime  
                  ├── momdip_mtx
                  ├── eigenvalues 
                  ├── eigenvectors
                  ├── transpose
                  ├── momdip_es_mtx               
                  └── source_file

Note that a copy of the source file has also been created into the ``data`` directory.

.. Hyperlink targets

.. _linalg.eig: https://numpy.org/doc/stable/reference/generated/numpy.linalg.eig.html
.. _transpose: https://numpy.org/doc/stable/reference/generated/numpy.transpose.html
.. _NumPy: https://numpy.org/
.. _matmul: https://numpy.org/doc/stable/reference/generated/numpy.matmul.html#numpy.matmul
.. _Q-CHEM: https://www.q-chem.com/