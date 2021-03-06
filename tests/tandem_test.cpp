#include "../tandem.h"
#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <fstream>
#include <iostream> // std::cout
#include <math.h>   // sqrt
#include <string>
#include <vector>
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE HartreeFockTest
#include <boost/test/included/unit_test.hpp>  // include this to get main()

using namespace std;

// The most simple test function possible, ground state of the
// non-interacting system of 3 particles in 4 orbitals. |0111>. This
// test the switch case in witch the difference betwen states is 0.
BOOST_AUTO_TEST_CASE(zero_differences) {
  // System details
  const int num_orbitals = 4;
  const int num_particles = 3;
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);

  // Wave function
  vector<double> wave_function(basis_length, 0);
  wave_function.at(0) = 1;

  // Expected output
  vector<double> expected_two_matrix = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                       0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
  Tandem alg(num_orbitals, num_particles, wave_function);
  vector<double> two_matrix = alg.run();
  assert(two_matrix.size() == expected_two_matrix.size());
  for (unsigned long int i = 0; i < two_matrix.size(); i++) {
    if (abs(expected_two_matrix.at(i) - two_matrix.at(i)) > 0.0001) {
      cout << "FAILED!" << endl;
      cout << "GOT: ";
      for (unsigned long int j = 0; j < two_matrix.size(); j++) {
        cout << two_matrix.at(j) << " ";
      }
      cout << endl;
      cout << "Expected: ";
      for (unsigned long int j = 0; j < two_matrix.size(); j++) {
        cout << expected_two_matrix.at(j) << " ";
      }
      cout << endl;
      BOOST_REQUIRE(false);
    }
  }
}

BOOST_AUTO_TEST_CASE(two_differences) {
  // System details
  const int num_orbitals = 4;
  const int num_particles = 3;
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);

  // Wave function
  vector<double> wave_function(basis_length, 0);
  wave_function.at(0) = 1;
  wave_function.at(1) = 1; // screw normalization

  // Expected output
  vector<double> expected_two_matrix = {2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
                                       1, 0, 0, 0, 1, 1, 0, 1, 0, 0};

  Tandem alg(num_orbitals, num_particles, wave_function);
  vector<double> two_matrix = alg.run();
  assert(two_matrix.size() == expected_two_matrix.size());
  for (unsigned long int i = 0; i < two_matrix.size(); i++) {
    if (abs(expected_two_matrix.at(i) - two_matrix.at(i)) > 0.0001) {
      cout << "FAILED!" << endl;
      cout << "GOT: ";
      for (unsigned long int j = 0; j < two_matrix.size(); j++) {
        cout << two_matrix.at(j) << " ";
      }
      cout << endl;
      cout << "Expected: ";
      for (unsigned long int j = 0; j < two_matrix.size(); j++) {
        cout << expected_two_matrix.at(j) << " ";
      }
      cout << endl;
      BOOST_REQUIRE(false);
    }
  }
}

BOOST_AUTO_TEST_CASE(four_differences) {
  // System details
  const int num_orbitals = 5;
  const int num_particles = 3;
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);

  // Wave function
  vector<double> wave_function(basis_length, 0);
  wave_function.at(0) = 1;
  wave_function.at(7) = 1; // screw normalization

  // Expected output
  vector<double> expected_two_matrix = {
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, // <- This one is
                                                               // the case 4
                                                               // relevant one!
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  Tandem alg(num_orbitals, num_particles, wave_function);
  vector<double> two_matrix = alg.run();
  assert(two_matrix.size() == expected_two_matrix.size());
  for (unsigned long int i = 0; i < two_matrix.size(); i++) {
    if (abs(expected_two_matrix.at(i) - two_matrix.at(i)) > 0.0001) {
      cout << "FAILED!" << endl;
      cout << "GOT: ";
      for (unsigned long int j = 0; j < two_matrix.size(); j++) {
        cout << two_matrix.at(j) << " ";
      }
      cout << endl;
      cout << "Expected: ";
      for (unsigned long int j = 0; j < two_matrix.size(); j++) {
        cout << expected_two_matrix.at(j) << " ";
      }
      cout << endl;
      BOOST_REQUIRE(false);
    }
  }
}
