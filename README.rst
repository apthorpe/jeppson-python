==============
jeppson-python
==============


Pipe network flow analysis toolkit


Description
===========


A library and set of applications replicating the software in *Steady Flow
Analysis of Pipe Networks: An Instructional Manual* (1974). *Reports.* Paper
300.  http://digitalcommons.usu.edu/water_rep/300 and
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


Design Information
==================

All programs make use of the InputLine class (in jeppson.input) for tokenizing
input lines and differentiating between blank, comment, and input data lines.
The original Fortran applications could not process blank lines or '#'-prefixed
comment lines; the Python versions are substantially more robust in input
processing and validation.

Rather than implementing the Darcy-Weisbach friction factor correlation given
in the original Fortran applications, the jeppson_ch2, jeppson_ch4,
jeppson_ch5, and jeppson_ch7 programs all make use of the friction_factor()
from the fluids.friction library. While this makes each program more of a
'black box', hiding the complexity of friction factor calculation in the call
to an external library, friction factors are calculated from empirical
correlations which themselves are a different sort of black box. It's not clear
anything substantial is lost in using an external library versus independently
coding the friction factor correlation. The ``fluids`` library contains a large
number of correlations which had been developed after the Jeppson texts had
been published. In practice, the friction_factor() routine produces results
very similar to the original Fortran implementation as should be expected.

The Pint physical unit library was used in each program to ensure dimensional
consistency and allow for simplified unit conversions. These applications show
the distinction between *physical computing* and *numerical computing*.
Whenever possible, variables representing physical quantities are stored as
Pint Quantity objects and physical equations are performed using Quantity
objects as well. This removes the need to mannually apply conversion factors in
calculation or display. These programs serve as technical demonstrators for
using the Pint library for physical computing and analysis.

The Python applications make extensive changes to the data models used in the
original Fortan applications. The most advanced data structure used by the
Fortran programs is the array; in the Python versions, state data is stored in
a complex data structure composed of lists, dicts, and sets. The complex data
storage allow for more direct, more obvious access to data which (ideally)
makes the underlying theory of each program much clearer. Additionally, since
Python data structures are flexible, this should result in reduced memory use.

In programs jeppson_ch5, jeppson_ch6a, and jeppson_ch7, the pygraphviz library
was used to generate flow topology diagrams for displaying results and
validating input.


Note
====


This project has been set up using PyScaffold 3.0.3 via

    putup -p jeppson -d "Pipe network flow analysis toolkit" -l mit jeppson-python

For details and usage information on PyScaffold see http://pyscaffold.org/.
