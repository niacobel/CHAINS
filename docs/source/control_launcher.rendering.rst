***********************
Rendering the templates
***********************

What is rendering?
==================

.. todo::

   COMING SOON

Jinja templates
===============

Rather than trying to explain how to create your templates in a vacuum, let's consider two basic examples to illustrate this process.

.. warning::

   All the Jinja templates must be placed inside the ``templates`` directory of ``CONTROL LAUNCHER``.

Parameters file template
------------------------

.. todo::

   COMING SOON

Job script template
-------------------

.. todo::

   COMING SOON

YAML configuration file
=======================

.. todo::

   COMING SOON

.. _control_rendering_fct:

Rendering functions
===================

.. todo::

   COMING SOON

General definition
------------------

.. todo::

   COMING SOON

The six arguments
-----------------

.. todo::

   COMING SOON

Using the Jinja templates
-------------------------

In order to use and render the Jinja templates, the rendering functions call another function, named ``jinja_render`` and defined in the ``control_renderer.py`` file. This function is the one that makes the link between the templates and their specific content:

.. autofunction:: control_renderer.jinja_render

This function receives three key arguments describing where to find the template and what is the corresponding content. It then returns the content of the rendered file in a variable (``output_text``), that will be printed in a file by the main script ``control_launcher.py``.

Simple function model
---------------------

.. todo::

   COMING SOON

Calling your rendering function
-------------------------------

.. todo::

   COMING SOON

Example of a rendering function
-------------------------------

.. todo::

   COMING SOON

Review of the template and configuration files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo::

   COMING SOON

Function definition
~~~~~~~~~~~~~~~~~~~

.. todo::

   COMING SOON

.. Hyperlink targets
