// The Energy test. The mechanism with which we find negative examples
// is by calculating the energy expectation value with respect to some
// solved hamiltonians. If the expectationvalue is below the known
// ground state energy then the 2-matrix cannot be N-Representable.

// This program takes the energy expectationvalue of num_examples
// generated 2-matrices with each of the num_hamiltonoans hamiltonians
// that are read from file.

#include "../tandem.h"
#include "test.h"
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp> // not sure if needed
#include <fstream>
#include <iostream> // std::cout
#include <math.h>   // sqrt
#include <random>
#include <string>
#include <vector>

using namespace std;

double abs(vector<double> vec) {
  double res = 0;
  for (int i = 0; i < vec.size(); i++) {
    res += vec.at(i) * vec.at(i);
  }
  return sqrt(res);
}

int main() {
  cout << "energy test... ";
  // Basic info
  const int num_orbitals = 10;      // Dependent on hamiltonians!
  const int num_particles = 4;      // Dependent on hamiltonians!
  const int num_hamiltonians = 100; // Dependent on hamiltonians!
  const string file_name_hamiltonians = "../../../data/H.dat";
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  const int num_examples = 10;
  // Initialize vectros
  vector<double> wave_function(basis_length);
  vector<vector<double>> hamiltonians(num_hamiltonians,
                                     vector<double>(pow(num_orbitals, 4), 0));
  // reading in hamitonians to test system (4 particles 10 orbitals)
  ifstream file_hamiltonians(file_name_hamiltonians, ios_base::in);
  // Check if file opended succesfully
  if (!file_hamiltonians) {
    cout << "ERROR: energy test: " << file_name_hamiltonians
         << " did not open succesfully" << endl;
    return -1;
  }

  int counter = 0;
  vector<double> energy(num_hamiltonians), interaction(num_hamiltonians);
  for (std::string line;
       std::getline(file_hamiltonians, line) and counter < num_examples;
       counter++) // read stream line by line
  {
    std::istringstream in(line); // make a stream for the line itsel
    in >> energy.at(counter) >> interaction.at(counter);
    int counter2 = 0;
    double element;
    while (in >> element) {
      hamiltonians.at(counter).at(counter2) = element;
      counter2++;
    }
    counter2 = 0;
  }
  // Define a random number generator and initialize it with a reproducible
  // seed.
  // (The seed is unsigned, otherwise the wrong overload may be selected
  // when using mt19937 as the base_generator_type.)
  boost::minstd_rand generator(42u);

  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::uniform_real<double> uni_dist(-1,1);
  boost::variate_generator<boost::minstd_rand&, boost::uniform_real<double> > uni(generator, uni_dist);
  // You can now retrieve random numbers from that distribution by means
  // of a STL Generator interface, i.e. calling the generator as a zero-
  // argument function
  for (int i = 0; i < num_examples; i++) {
    // generate uniformly random wave functions
    for (int j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = uni();
    }
    // Normalize wave function
    double total_sum = abs(wave_function.at(i));
    for (int j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = wave_function.at(j) / total_sum;
    }
    // Initialize and run algorithm
    Tandem alg(num_orbitals, num_particles, wave_function);
    vector<double> two_matrix = alg.run();
    Test test(alg);
    for (size_t j = 0; j < num_hamiltonians; j++) {
      double E = test.energy(hamiltonians[j]);
      if (energy.at(j) > E) {
        cout << "ERROR: energy test FAILED!" << endl;
      }
    }
  }
  cout << " PASSED!" << endl;
  return 0;
}
