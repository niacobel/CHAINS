***************************
Determining the transitions
***************************

To determine the transitions that will be covered by the control procedure, ``CONTROL LAUNCHER`` calls a function named the **transition function** and defined in the ``transition_fcts.py`` file. This function creates the initial and target states file (also sometimes called "population files") and defines the different combinations of initial and target states (a.k.a the transitions). The only arguments this function needs are the ``system`` dictionary and the path towards the ``data`` directory, where the states file will be created. If you do not know what are those two variables, please refer to the :doc:`control_launcher.modelling` specific documentation.

General definition of the transition functions
==============================================

All the transition functions are defined in the ``transition_fcts.py`` file and obey some restrictions, in order to be callable by ``CONTROL LAUNCHER``:

- They must only take two arguments: the ``system`` dictionary and the path towards the ``data`` directory.
- They must return a list of dictionaries (``transitions_list``) containing one mandatory key, ``label``, which is the name of the considered transitions, used for naming the job directories. You are of course likely to have additional keys, such as the name of the created files that may be needed by your :ref:`rendering function <control_rendering_fct>`, but this is the only one explicitly needed by the main script itself.

If a problem arises when determining the transitions, a ``ControlError`` exception should be raised with a proper error message (see :ref:`how to handle errors <control_errors>` for more details).

Choosing a transition function
==============================

The transition function that will be called by ``CONTROL LAUNCHER`` is the one associated with the ``transition_function`` YAML key defined in the :ref:`clusters configuration file <control_clusters_file>`:

.. code-block:: yaml

   mycluster:
      profiles:
         myprofile1:
            transition_function: name-of-transition-function
         myprofile2:
            transition_function: name-of-transition-function

where ``mycluster`` corresponds to the name of your cluster (given as a :ref:`command line argument <control_arguments>`) while ``myprofile1`` and ``myprofile2`` are the names of the profiles you want to run (such as ``chains_qoctra``). This way, a different transition function can be assigned to each profile.

proj_ground_to_triplet transition function
==========================================

.. autofunction:: transition_fcts.proj_ground_to_triplet

The transitions defined by this function are the ones that transfer the whole population from the ground state to each of the triplet states, individually. To identify all the triplet states of the molecule, this function uses the information contained in the additional ``states_list`` key defined in ``system`` by the parsing function.

This function thus creates one initial state file, named ``ground_1``, where the whole population lies in the ground state, and a number of target state files equal to the number of triplets in the molecule, where the whole population lies in one of the triplets. Those target states files are defined as projectors (used by the ``OPM`` operating mode of QOCT-RA) and are named ``projector<triplet_label>_1``, where *<triplet_label>* is the value of the ``label`` key of the triplet state in the ``states_list`` dictionary.  All the files are created as density matrices.

.. note::

   Another, dummy file is created by this function, named ``final_1``. This file does not contain any relevant information but is still needed by QOCT-RA in order to run.

As an example, let's use the |H2O| case presented in the :doc:`control_launcher.modelling` specific documentation. As a reminder, this is the states list of this molecule:

.. |H2O| replace:: H\ :sub:`2`\ O

.. code-block:: text

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

At the end of this process, the output directory structure is then

.. code-block:: text

    out_dir/ 
      └── source/
            └── data/
                  ├── all the data files created during the system modelling step
                  ├── ground_1
                  ├── projectorT1_1
                  ├── projectorT2_1
                  ├── projectorT3_1
                  ├── projectorT4_1 
                  └── final_1

The content of ``transitions_list`` in this case is

.. code-block:: python

   transitions_list = [
      {"label" : T1, "init_file" : ground_, "target_file" : projectorT1_},
      {"label" : T2, "init_file" : ground_, "target_file" : projectorT2_},
      {"label" : T3, "init_file" : ground_, "target_file" : projectorT3_},
      {"label" : T4, "init_file" : ground_, "target_file" : projectorT4_}
      ]

And here is the content of the created files:

.. code-block:: text
   :caption: ground_1

   ( 1.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )

.. code-block:: text
   :caption: projectorT1_1 (T1 is the 2nd state)

   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 1.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )

.. code-block:: text
   :caption: projectorT2_1 (T2 is the 4th state)

   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 1.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )

.. code-block:: text
   :caption: projectorT3_1 (T3 is the 5th state)

   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 1.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )

.. code-block:: text
   :caption: projectorT4_1 (T4 is the 7th state)

   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 1.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )

.. code-block:: text
   :caption: final_1

   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )
   ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 ) ( 0.00 , 0.00 )