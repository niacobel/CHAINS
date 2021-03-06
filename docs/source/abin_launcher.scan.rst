**************************
Scanning the geometry file
**************************

While ``ABIN LAUNCHER`` scans the geometry file, it is looking for information about the chemical formula and the atomic coordinates of our molecule. The way it does that is by calling a function defined in the ``geom_scan.py`` file, called a **scanning function**. The only argument this function needs is the content of the geometry file, stored as a list in the ``mol_content`` variable (each element of the list corresponds to a line of the file). 

Note that at this time, only the XYZ format is supported for geometry files and so there is only a single scanning function, called ``xyz_scan``. The XYZ format is rather well known and you should be able to convert your geometry files to this format quite easily. However, if the need arises, consult the :ref:`other_formats` subsection for more details on how to handle another format.

XYZ format
==========

Structure of XYZ geometry files
-------------------------------

Let's take a geometry file for the |SiH4| molecule as an example:

.. code-block:: text

   5
   my comment
   Si        -0.31438        1.89081        0.00000
   H          1.10562        1.89081        0.00000
   H         -0.78771        2.04313       -1.33010
   H         -0.78772        2.96655        0.79696
   H         -0.78772        0.66276        0.53314

The structure of this file is as follows:

- The first line of the file is the number of atoms in the molecule.
- The second line is a comment line where you can write anything you want. 
- The rest of the file is the list of the atomic coordinates. They are written in the form *<atom type> <x coordinate> <y coordinate> <z coordinate>*. 

.. |SiH4| replace:: SiH\ :sub:`4`\ 

Function definition
-------------------

.. autofunction:: geom_scan.xyz_scan

Defined in the ``geom_scan.py`` file, the ``xyz_scan`` function receives the content of the geometry file as an argument (a list variable) and returns a ``file_data`` dictionary with two keys: ``chemical_formula`` and ``atomic_coordinates``: 

- The first key, ``chemical_formula``, contains a dictionary where each *key:value* pair corresponds to an atom type and the number of atoms of that type in the molecule. For example, the value for the ``chemical_formula`` key for |H3PO4| is *{'H':3, 'P':1, 'O':4}*.
- The second key, ``atomic_coordinates`` is simply a list of the atomic coordinates lines of the XYZ file.

By making use of :ref:`regular expressions <regex>`, the function detects all the lines that have the *<atom type> <x coordinate> <y coordinate> <z coordinate>* pattern and stores them as a list into ``atomic_coordinates``. For ``chemical_formula``, the function counts the number of occurrences of each atom type as it reads it at the beginning of each coordinates line.

.. |H3PO4| replace:: H\ :sub:`3`\ PO\ :sub:`4`\ 


.. _other_formats:

Handling other formats
======================

If you want to add support for another format, you will need two things:

- Add a scanning function for that format to ``geom_scan.py``
- Tell ``ABIN LAUNCHER`` to use that new function, most likely by adding a new command line argument

Format as a new command line argument
-------------------------------------

If you want to add support for another format, it is probably a good idea to add a new command line argument that will be given when executing ``ABIN LAUNCHER``. Everything has already been prepared for that eventuality. As such, you only need to change two lines in ``abin_launcher.py``:

Just uncomment this line (remove the # at the beginning of the line)

.. code-block:: python

   required.add_argument('-f', '--format', type=str, help="Format of the geometry files that need to be read.", required=True)

then replace this line

.. code-block:: python

   mol_fmt = "xyz"

by

.. code-block:: python

   mol_fmt = args.format

and it's done! 

Now you will be able to specify the format of your geometry files by using the :guilabel:`-f / \\--format` command line argument when executing ``abin_launcher.py``.

.. warning::

   The extension of the geometry files must be the same as the value given in the :guilabel:`-f / \\--format` argument.

Defining a new scanning function
--------------------------------

All the scanning functions must be defined in the ``geom_scan.py`` file and need to obey some restrictions in order to be callable by ``ABIN LAUNCHER``:

- They must be named *fmt_scan*, where *fmt* is the format of the geometry file as it will be given in the command line.
- They must only take one argument: a list containing the lines of the geometry file (``mol_content``).
- They must return a dictionary (``file_data``), following the pattern :code:`{ 'chemical_formula' : { }, 'atomic_coordinates' : [ ] }`. You can have additional keys if you want, but those two are mandatory.
  
If a problem arises when scanning the molecule file, an ``AbinError`` exception should be raised with a proper error message (see :ref:`how to handle errors <abin_errors>` for more details).