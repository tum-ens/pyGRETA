.. renewable-timeseries documentation master file, created by
   sphinx-quickstart on Fri Jul 26 11:33:25 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Renewable-timeseries: GIS-based model for renewable energy potential prediction and timeseries genertion
========================================================================================================

:Author: Kais Siala <kais.siala@tum.de>, Houssame Houmy <houmyh@gmail.com>
:Organization: `Chair of Renewable and Sustainable Energy Systems`_,
               Technical University of Munich,
:Version: |version|
:Date: |today|
:Copyright:
  The model code is licensed under the `GNU General Public License 3.0
  <http://www.gnu.org/licenses/gpl-3.0>`_.
  This documentation is licensed under a `Creative Commons Attribution 4.0 
  International <http://creativecommons.org/licenses/by/4.0/>`_ license.

Contents
--------

User's manual
^^^^^^^^^^^^^

We can link here to the user manual for the config files maybe

.. the following section contains the links to the other parts of the documentation, the 'maxdepth' component define how many sub sections should be displayed

.. toctree::
   :maxdepth: 1
	

Mathematical documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Continue here if you want to understand the theoretical conception of the model
generator, the logic behind the equations and the structure of the features.

.. toctree::
   :maxdepth: 1

   theoretical

Technical documentation
^^^^^^^^^^^^^^^^^^^^^^^

Continue here if you want to understand in detail the model generator
implementation.

.. toctree::
   :maxdepth: 1
   
   implementation


Features
--------
* To be completed 


Dependencies
------------

* `Python`_ versions 2.7 or 3.x are both supported.
* `pyomo`_ for model equations and as the interface to optimisation solvers.
* `pandas`_ for input and result data handling, report generation ...
* Any solver supported by pyomo; suggestion: `GLPK`_,gurobi
* To be completed...

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
   
.. _glpk: https://www.gnu.org/software/glpk/
.. _Chair of Renewable and Sustainable Energy Systems: http://www.ens.ei.tum.de/
.. _matplotlib: http://matplotlib.org
.. _pandas: http://pandas.pydata.org
.. _pyomo: http://www.pyomo.org
.. _python: https://www.python.org/
.. _readme.md: https://github.com/tum-ens/urbs/blob/master/README.md#installation
.. _urbs: https://github.com/tum-ens/urbs
