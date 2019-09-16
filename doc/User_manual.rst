User manual
===========

Installation
------------

.. NOTE:: We assume that you are familiar with *git* and *conda*.

First, clone the git repository in a directory of your choice using a Command Prompt window::

	$ ~\directory-of-my-choice> git clone https://github.com/tum-ens/renewable-timeseries.git

We recommend using conda and installing the environment from the file ``PyMTS37.yml`` that you can find in the repository. In the Command Prompt window, type::

	$ cd renewable-timeseries\env\
	$ conda env create -f PyMTS37.yml

Then activate the environment::

	$ conda activate PyMTS37

In the folder ``code``, you will find multiple files:

+-----------------------+---------------------------------------------------------------------------+
| File                  | Description                                                               |
+=======================+===========================================================================+
| config.py             | used for configuration, see below.                                        |
+-----------------------+---------------------------------------------------------------------------+
| Master.py             | main file, which will be run later using ``python Master.py``.            |
+-----------------------+---------------------------------------------------------------------------+
| data_functions.py     | contains helping functions for data formatting, filtering, reshaping, etc.|
+-----------------------+---------------------------------------------------------------------------+
| model_functions.py    | contains helping functions for the modeling.                              |
+-----------------------+---------------------------------------------------------------------------+
| util.py               | contains minor helping functions and the necessary python libraries to be |
|                       | imported.                                                                 |
+-----------------------+---------------------------------------------------------------------------+

config.py                                                                                           
---------
This file contains the user preferences, the links to the input files, and the paths where the outputs should be saved.
The paths are initialized in a way that follows a particular folder hierarchy. However, you can change the hierarchy as you wish.

.. toctree::
   :maxdepth: 3
   
   source/config


Recommended input sources
-------------------------

Expected outputs
----------------

Recommended workflow
--------------------
The script is designed to be modular, yet there is a recommended work flow to follow for your first run...
