#include "../tandem.h"
#include "../cholesky.h"
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp> // not sure if needed
#include <boost/math/special_functions/binomial.hpp>
#include <boost/filesystem.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream> // std::cout
#include <math.h>   // sqrt
#include <random>
#include <string>
#include <vector>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <functional>
#include <ctime>
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE BenchmarkTest
#include <boost/test/included/unit_test.hpp>  // include this to get main()


double abs(vector<double> vec) {
  double res = 0;
  for (size_t i = 0; i < vec.size(); i++) {
    res += vec.at(i) * vec.at(i);
  }
  return sqrt(res);
}

BOOST_AUTO_TEST_CASE( hartree_fock_test ) {
  clock_t tStart = clock();
    // System details
  const int num_orbitals = 10;
  const int num_particles = 4;
  const int num_examples = 100;
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  // Define a random number generator and initialize it with a reproducible
  // seed.
  // (The seed is unsigned, otherwise the wrong overload may be selected
  // when using mt19937 as the base_generator_type.)
  boost::minstd_rand generator(42);

  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::uniform_real<double> uni_dist(-1,1);
  boost::variate_generator<boost::minstd_rand&, boost::uniform_real<double> > uni(generator, uni_dist);
  // You can now retrieve random numbers from that distribution by means
  // of a STL Generator interface, i.e. calling the generator as a zero-
  // argument function
  vector<double> wave_function(basis_length);
  for (size_t i = 0; i < num_examples; i++) {
    for (size_t j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = uni();
    }
    double total_sum = abs(wave_function);
    for (size_t j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = wave_function.at(j) / total_sum;
    }
    assert(abs(abs(wave_function) - 1) < 0.001);
    Tandem alg(num_orbitals, num_particles, wave_function);
    vector<double> two_matrix = alg.run();
  }
 BOOST_REQUIRE((double)(clock() - tStart)/CLOCKS_PER_SEC < 3);
}
