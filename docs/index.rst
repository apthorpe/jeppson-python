==============
jeppson-python
==============

Pipe network flow analysis toolkit based on the work of Roland W. Jeppson.

``python-jeppson`` is a library and set of applications replicating the
software in *Steady Flow Analysis of Pipe Networks: An Instructional Manual*
(1974). *Reports.* Paper 300.  http://digitalcommons.usu.edu/water_rep/300 and
*Analysis of Flow in Pipe Networks* (1976). Ann Arbor Science Publishers, Inc.
http://www.worldcat.org/title/analysis-of-flow-in-pipe-networks/oclc/927534147
by Roland W. Jeppson.

Six command line applications are included, providing the functionality of the
original Fortran applications described in the Jeppson texts:

* ``jeppson_ch2`` - Frictional head loss calculator
* ``jeppson_ch4`` - Incompressible flow calculator
* ``jeppson_ch5`` - Linear method pipe network flow solver
* ``jeppson_ch6a`` - Newton-Raphson solver which determines junction pressures
* ``jeppson_ch6b`` - Newton-Raphson solver which generates corrective loop flows
* ``jeppson_ch7`` - Hardy Cross method pipe network flow solver

Each program takes the same input file structure as its Fortran equivalent and
generates the same results, though the output format may differ substantially
from the original Fortran application.

The original Fortran applications were recovered and modernized in a separate
project, located at https://bitbucket.org/apthorpe/jeppson_pipeflow A brief
description of the recovery process and a whitepaper summarizing the insights
gained from the recovery and modernization work can be found at
https://www.linkedin.com/pulse/case-study-revitalizing-legacy-engineering-jeppson-pipe-bob-apthorpe/

Consider this project to be thematically related, demonstrating the
implementation of these programs in Python.

Contents
========

.. toctree::
   :maxdepth: 2

   License <license>
   Authors <authors>
   Changelog <changelog>
   Module Reference <api/modules>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
