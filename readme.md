Welcome TANDEM - TrAce N DEnsity Matrices

This module is all about generating examples to train on. Many
different ways to generate. Tandem is also probably the best
maintained package between (annan, scripts and data) as it is most
likley not to be completley irrelevant.

The python bindings, documented in PyTandem.cpp are the most usefull
and expose allmost all of the packages functionality. The rest is
legacy code and tests.

To make it easier to maintain and make changes this module has a
.travis.yml file that facilitates automatic unittesting when commiting
to a travis enable git repository.

# Installation:

Create a python2 virtual enviroment and activate it then do:

  - pip install -r requirements.txt
  - export PYTHONPATH=$PWD
  - sudo apt-get install libboost-all-dev libgsl0-dev
  - make
  - make -C tests

To compile and run tests. After this the functions ported in
PyTandem.cpp will be importable from PyTandem.

This should work most of the time on ubuntu. Depending on you
installation some more packages might have to be installed via the
package manage. On other distros the apt-get command will fail.

If you are using windows: Create an ubuntu virtual machine or this
code base will be most probably be a living hell to use.




