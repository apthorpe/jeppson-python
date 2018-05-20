==============
jeppson-python
==============


Pipe network flow analysis toolkit based on the work of Roland W. Jeppson.


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

The original Fortran applications were recovered and modernized in a separate
project, located at https://bitbucket.org/apthorpe/jeppson_pipeflow A brief
description of the recovery process and a whitepaper summarizing the insights
gained from the recovery and modernization work can be found at
https://www.linkedin.com/pulse/case-study-revitalizing-legacy-engineering-jeppson-pipe-bob-apthorpe/

Consider this project to be thematically related, demonstrating the
implementation of these programs in Python.


Disclaimer
==========


This software is provided as an educational resource and no warranty is made to
its accuracy or suitability for engineering use, especially in any use
involving protection of human life and environmental quality. *Caveat
utilitor!*


Project Goals
=============

The primary goals of this project were to write Python equivalents to the six
command-line interface (CLI) applications described in *Steady Flow Analysis of
Pipe Networks: An Instructional Manual* (1974) and *Analysis of Flow in Pipe
Networks* (1976); see above. The programs should be able to process the input
files used by the original Fortran applicaitons and should produce (at minimum)
the same output (content, precision, accuracy). The Python applications should
be written in standard Python idiom (i.e. should not be a transliteration of
the original Fortran into Python), and should demonstrate the use of diagnostic
logging, standard CLI argument and option processing, error- and exception
handling, and user-centric error messages. If possible, pipe network topology
and pipe flows should be presented visually using GraphViz -
https://www.graphviz.org/

The Python applications are not required to follow the implementation methods
used in the original Fortran; it is expected that full advantage will be taken
of Python's data structures and coding constructs as well as standard and
third-party libraries.

Similarly, the output produced by the Python applications may vary
substantially from the original Fortran applications provided the same content
may be found in the output of the new applications.

A major goal is to demonstrate the use of the Pint unit conversion and physical
quantity library. Since poor unit conversion and specification have led to a
number of dramatic and expensive failures of mission- and safety-critical
software, the project serves to demonstrate the use and limitations of
unit-aware physical computing.

Finally, a goal in rewriting the network flow solvers was to iteratively refine
and 'standardize' applications to find common elements, structures, or methods
which could be extracted into resuable components. It was not clear at the
outset if a flow solver object model would suggest itself in the course of
writing the Python applications so this project was partly intended to find
code duplicated between applications which could be refactored into independent
objects or libraries. Rather than proposing an object model early in the
project which may or may not serve the application, the applications were
intentionally developed in a 'naive' manner, hoping similarities between them
would suggest evolutionary refactoring and componentization. This allowed focus
to be put primarily on generating correct results, secondarily on
object-oriented design. This focus allowed for rapid development and provided
test cases to detect any errors which appeared during refactoring.


Design Information
==================


The underlying theory of operation and input file format description for each
of the applications can be found in the appropriate chapter of the Jeppson
documents, for example http://digitalcommons.usu.edu/water_rep/300

All command-line applications support the ``-v`` and ``-vv`` flags to increase
the verbosity of diagnostic information shown. This is implemented via Python's
standard logging module. Logging is extended to the underlying libraries
(jeppson.input).

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

An object-oriented approach toward the full application was not taken since
there was little chance of code reuse. A data model was defined specifically
for each of the network flow solver programs (jeppson_ch5, jeppson_ch6a,
jeppson_ch6b, and jeppson_ch7). The data models have similar sections and
overal structure, but they are not interchangeable. Example data models for
each network flow solver can be found in the ``userdoc`` directory.

In programs jeppson_ch5 and jeppson_ch6a, the pygraphviz library was used to
generate flow topology diagrams for displaying results and validating input. At
present the diagram layout and quantitative elements need tuning to improve
presentation, however diagrams are complete and accurate with respect to
topology, pressure and flow display, inflows, outflows, and flow direction.

Unit testing is applied principally to object-oriented components, mainly
the Pipe class in jeppson.pipe and the InputLine class in jeppson.input.
Integral testing was used to manually compare the Python applications with the
original Fortran applications. Had the applications been designed as objects,
unit testing would have been a more reasonable choice since it can easily be
automated. The choice of application architecture makes the individual programs
rather difficult to test; a more modular or object-oriented design would
simplify testing but would also complicate implementation. In this case, the
decision was to go with a simpler application architecture and trade ease of
implementation for ease of testing. This is reasonable in a prototype or
demonstrator application such as this; it may not be appropriate for other
application roles and use cases.


Possible Future Work
====================


Add GraphViz support for all network solvers
--------------------------------------------

Adding GraphViz support for the ``jeppson_ch5`` and ``jeppson_ch6a``
applications was straightforward since the data model for these solvers
included lists of pipes and the junctions which connect them.
GraphViz support for ``jeppson_ch6b`` and ``jeppson_ch7`` programs is
complicated since junctions are not explicitly enumerated; instead, a list of
pipes and flow loops comprised of pipes are provided. Generating and
enumerating a list of connecting junctions from loop data is a fairly involved
process.

Data serialization
------------------

Serializing the case_dom data model used in the network flow solvers in a
format such as YAML or JSON would simplify post-processing the code results.

Structured input
----------------

Converting from free-form text input to a serialized input format such as JSON
or YAML would allow the code to be driven with a different user interface (e.g.
web, desktop GUI)

The case_dom data structure may be more useful if converted to several
independent Pandas data frames, then joined or queried in order to simplify
data access. This may be useful both for matrix and vector construction while
solving for network flows or for post-processing, analysis, and visualization.

Improved data visualization
---------------------------

The network flow solvers produce a conservative directed graph of volumetric
flow, ideal for representation in a Sankey plot; see
https://www.sciencedirect.com/science/article/pii/S0921344917301167?via%3Dihub


Documentation
=============

The code has been developed with the intent of using Sphinx as the project
documentation processor; http://www.sphinx-doc.org/en/master/ All files,
classes, and functions should have the appropriate docstring present. Functions
will additionally describe required and optional arguments, return values, and
exceptions raised, if any.

Documentation of the code theory or input format are available in the original
Jeppson documents and are not repeated here. This is considered reasonable
since this project is part of a larger whole based on Jeppson's original
texts.


Testing
=======

Testing is discussed near the end of the *Design Information* section. Unit
testing has been used extensively on component classes (InputLine and Pipe
classes). The command-line applications are primarily tested via integral
testing to ensure existing input files may be read and results of the original
and Python applications are comparable and reasonably close (within 5-10%).
Identical numerical results are not expected due to differences in precision
and calculational method. Note the discussion of design trade-offs and ease of
testing in the *Design Information* section.

The network flow solvers (``jeppson_ch5``, ``jeppson_ch6a``, ``jeppson_ch6b``,
and ``jeppson_ch7``) are all tested via integral tests using the pytest
framework, so it is possible to automate integral testing. It is not as simple
as unit testing (relying on several fixtures) and the fidelity is rather coarse
but it can be done.

Additionally, flake8 compliance is incorporated in the test suite to enforce a
reasonable level of stylistic quality.


Development Note
================


This project has been set up using PyScaffold 3.0.3 via

    putup -p jeppson -d "Pipe network flow analysis toolkit" -l mit jeppson-python

For details and usage information on PyScaffold see http://pyscaffold.org/.
