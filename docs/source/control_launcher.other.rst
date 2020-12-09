*********************
Other important files
*********************

.. _control_clusters_file:

Clusters configuration file
===========================

.. todo::

   COMING SOON

.. _control_errors:

Errors handling
===============

When adding a :ref:`rendering function <control_rendering_fct>` or another custom function to ``CONTROL LAUNCHER``, having a way to handle errors is definitely useful. In ``CONTROL LAUNCHER``, this is managed by the ``control_errors.py`` file. It is somewhat basic but should be enough to cover your needs.

Custom exception
----------------

A custom exception class has been created to handle errors specific to ``CONTROL LAUNCHER``, in the ``control_errors.py`` file:

.. autoclass:: control_errors.ControlError 

Feel free to raise it when you want to prevent predictable errors from happening (missing file, incorrect value, etc.) by simply using

.. code-block:: python

   raise control_errors.ControlError ("my message here")

Those raised exceptions wil be caught by ``CONTROL LAUNCHER``, which will then either abort the execution or skip the incriminated configuration file, depending on where the error occurred.

Checking the existence of files and directories
-----------------------------------------------

In order to easily check if specific files or directories exist, the ``check_abspath`` function has been defined in the ``control_errors.py`` file:

.. autofunction:: control_errors.check_abspath

As an example, let's say we want to check if our clusters configuration file is still there, we can use the code:

.. code-block:: python

  clusters_file = control_errors.check_abspath(os.path.join(code_dir,"clusters.yml"),"YAML clusters configuration file","file")

This will check if there is a file named ``clusters.yml`` in ``CONTROL LAUNCHER``'s directory (``code_dir``) and if it is indeed a file (and not a directory for example). 

- If there is, it will return the absolute path towards that file (useful for referencing that file later in the script, no matter where the current directory is). 
- Otherwise, it will raise an exception and specify the context as "YAML clusters configuration file" for easy tracking, which will result in an error message of the form:

.. code-block:: console

   Something went wrong when checking the path  ~/CHAINS/control_launcher/clusters.yml
   Context:  YAML clusters configuration file
   ERROR: ~/CHAINS/control_launcher/clusters.yml does not seem to exist.

.. Hyperlink targets
