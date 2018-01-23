 TANDEM - TrAce N DEnsity Matrices

CLI frontend to the RA (Reduction Algorithm) class that reduces
N-particle wave functions to 2-particle density matrices.

USAGE: ./tandem -p number_of_particles -e number_of_examples -f file_to_save_to -d distribution

Implemented distributions are uniform, one_over_x-p and exponential-T
where p and T can be any number specifying 1/x**p and exp(-x/T) respectivly.

number of orbitals need to be known
at compile time. It can be changed by editing TANDEM_cli.cpp:37. Then
compile again with "make cli".

The library is dependent on BOOST and possibly it needs c++11.
