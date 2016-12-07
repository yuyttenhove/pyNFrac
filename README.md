# pyNFrac

A Python library wrapper around Sven De Rijcke's neutral fraction tables.

This library heavily depends on Boost Python and NumPy, and cannot compile
without these dependencies.

To compile, get a fresh copy of this repository using `git`, set up the build
environment using `cmake PATH_TO_MAIN_FOLDER`, and compile using `make`. The
generated `.so` library can be directly imported into Python (from within the
folder containing it, or from anywhere if you add this folder to your PATH).

Unit tests can be run using `make check`.
