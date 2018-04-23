#TANDEM - TrAce N DEnsity Matrices

Reduces N-particle wave functions to 2-particle density matrices.

USAGE: tandem -p number_of_particles -e number_of_examples -f file_to_save_to -d distribution

Implemented distributions are uniform, one_over_x-p and exponential-T
where p and T can be any number specifying 1/x^p and exp(-x/T) respectivly.

The library is dependent on BOOST and possibly it needs c++11.
