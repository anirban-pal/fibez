# fibez
This is Fibez, a tool to simulate the energetics and dynamics of fibrous materials.

FIBEZ - Fibers as cubic Beziers
===============================

This is Fibez, a tool to simulate the energetics and dynamics of fibrous materials.

Fibez is free software, you can redistribute it and/or modify it under
the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This package is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Availability
============

The current stable version of Fibez is available from Github:

Installation
============

The Fibez package has the following dependencies, with preferred versions indicated.

NLopt 2.7.1 (https://nlopt.readthedocs.io/en/latest/)
GSL 2.8 (https://www.gnu.org/software/gsl/)
Cuba 4.2.2 (https://feynarts.de/cuba/)

Please ensure that the above libraries are installed and the respective directories are provided in the Makefile. Installation can be performed by running "make" in the downloaded folder, which will create an executable "fibez".

Running the program
===================

The program "fibez" can be run from the command line with a single argument corresponding to a data file. The data file contains the topology of the fibrous network in terms of its control points. See sample data file (data.test) for details.

Other files can be modifed based on simulation needs. The file 'potential.h' may be modified if a different 2-body cohesive potential is needed. Currently, an LJ type potential is used. The file 'headers.h' may also be modified to change the simulation parameters.

To run the program, simply type
./fibez data.test

Program output
==============

The program will produce 3 different outputs. The first output (log.fibez) is a time series record of potential and kinetic energies and momenta during time integration and minimization. The other 2 outputs (Output_cps.lammpstrj, Output.lammpstrj) are LAMMPS format dump files containing control point and fiber trajectories. One may use a software tool like VMD or Ovito to visualize these files.

Reporting Bugs
==============

Please report bugs and issues via GitHub.
