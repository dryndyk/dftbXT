******************
Change Log DFTB+XT
******************

Notable project changes from release 1.01 (2018-01-18).
Changes in base DFTB+ are below.


Unreleased
==========

Added
-----


Changed
-------


Fixed
-----


1.02 (2019-01-18)
=================

Added
-----

- New git repository at GitHub github.com/tranas-open/dftbXT
  and re-fork from github.com/dftbplus/dftbplus/master.

- Integration with CP2K+XT.

- New web site tranas.org/opensuite/dftb+xt.html.


Changed
-------

- Write/Read options for contact Green functions and self-energies.


Fixed
-----


1.01 (2018-01-18)
=================

Added
-----

- Model Hamiltonians for transport calculations.
  We introduced the possibility to read model Hamiltonians from external data files and use it with
  or without a geometry structure. 

- Elastic dephasing.
  Two models of elastic dephasing ("Büttiker probe" and "vibronic dephasing") can be used now 
  to include the dephasing and dissipation beyond the coherent Green function method. 

- Application to STM spectroscopy.
  We added new options to simplify and make faster the calculation of currents for systems with 
  changeable geometry (like the STM setup). 


Changed
-------


Fixed
-----


****************
Change Log DFTB+
****************

Notable project changes since release 1.3.1 (2017-02-22).


Unreleased
==========

Added
-----

- Use of the ELSI library.
  
- Ability to perform ground state MD with excitation energies.

- Caching for transition charges in excited state.

- DFTB+ can be compiled as a library and accessed via high level API (version
  number is in the file api/mm/API_VERSION below the main directory).

- Onsite corrected hamiltonian for ground and Casidat excited state energies.

Changed
-------


Fixed
-----


18.2
====

Added
-----

- Option for removing translational and rotational degrees of freedom in modes.

- H5 correction for hydrogen bonds.


Changed
-------

- Updated parser version to 6.

- Syntax for H5 and DampedHX corrections for hydrogen bonds unified.


Fixed
-----

- Compilation when socket interface disabled.

- Stress tensor evaluation for 3rd order DFTB.

- Tollerance keyword typo.

- Corrected erroneous Lennard-Jones-dispersion for periodic cases (broken since
  release 1.3)

- Forces/stresses for dual spin orbit.


18.1 (2018-03-02)
=================

Added
-----

- MPI-parallelism.

- Various user settings for MPI-parallelism.

- Improved thread-parallelism.

- LBGFS geometry driver.

- Evaluation of electrostatic potentials at specified points in space.

- Blurred external charges for periodic systems.

- Option to read/write restart charges as ASCII text.

- Timer for collecting timings and printing them at program end.

- Tolerance of Ewald summation can be set in user input.

- Shutdown possibility when using socket driver.

- Header for code prints release / git commit version information.

- Warning when downloading license incompatible external components.

- Tool straingen for distorting gen-files.


Changed
-------

- Using allocatables instead of pointers where possible.

- Change to use the Fypp-preprocessor.

- Excited state (non-force) properties for multiple excitations.

- Broyden-mixer does not use file I/O.

- Source code documentation is Ford-compatible.

- Various refactorings to improve on modularity and code clarity.


Fixed
-----

- Keyword Atoms in modes_in.hsd consider only the first specified entry.

- Excited window selection in Cassida time-dependent calculation.

- Formatting of eigenvalues and fillings in detailed.out and band.out

- iPI socket interface with cluster geometries fixed (protocol contains
  redundant lattice information in these cases).


17.1 (2017-06-16)
=================

Added
-----

- Add dptools toolkit.


Changed
-------

- Convert to LGPL 3 license.

- Restructure source tree.

- Streamline autotest suite and build system.


Fixed
-----

- Skip irrelevant tests that give false positives for particular compilation
  modes.

- Make geometry writing in gen and xyz files consistent.
