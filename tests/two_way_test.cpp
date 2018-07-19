// two_way_test.cpp - hartree fock state test.
// For any orthonormal matrix D one can create a change of basis such that
// a^d_1*a^d_2*...*a^d_N = D^k(1)_1*D^k(2)_2*...*D^k(N)_N
// a^d_k(1)*a^d_k(2)*...*a^d_k(N)
// then it turns out the one matrix is D*D^T.

// The 1-matrix can also be calculated from the 2-matrix. This is what
// this program does. It reads wave_function and one_matrix generated
// by the scheme described above from files HF.wf (wave function) and
// one_mat.m for the 1-matrix. Calculates the 2-matrix from the HF.wf
// wave function then calculates the 1-matrix from that and then
// compares the calculated 1-matrix with the 1-matrix on file
// one_mat.m

#include "../tandem.h"
#include "test.h"
#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <fstream>
#include <iostream> // std::cout
#include <math.h>   // sqrt
#include <random>
#include <string>
#include <vector>
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE TwoWayTest
#include <boost/test/included/unit_test.hpp>  // include this to get main()

using namespace std;

typedef boost::numeric::ublas::matrix<double> boost_matrix;

double abs(vector<double> vec) {
  double res = 0;
  for (int i = 0; i < vec.size(); i++) {
    res += vec.at(i) * vec.at(i);
  }
  return sqrt(res);
}

BOOST_AUTO_TEST_CASE( hartree_fock_test )
{
  // Basic system info
  const int num_orbitals = 10;
  const int num_particles = 4;
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  const int num_examples = 1;
  vector<vector<double>> wave_function(num_examples,
                                      vector<double>(basis_length, 0));
  const string file_name_wf = "HF.wf";
  const string file_name_one_matrix = "one_mat.m";
  // Open files
  ifstream file_wf(file_name_wf, ios_base::in);
  ifstream file_one_matrix(file_name_one_matrix, ios_base::in);
  // read wave function
  unsigned int counter = 0;
  for (std::string line; std::getline(file_wf, line) and counter < num_examples;
       counter++) // read stream line by line
  {
    std::istringstream in(line); // make a stream for the line itsel
    int counter2 = 0;
    double wf;
    while (in >> wf) {
      wave_function.at(counter).at(counter2) = wf;
      counter2++;
    }
    counter2 = 0;
  }
  // read 1-matrix
  vector<vector<double>> one_matrix_from_file(num_orbitals,
                                             vector<double>(num_orbitals, 0));
  counter = 0;
  for (std::string line;
       std::getline(file_one_matrix, line) and counter < num_orbitals;
       counter++) // read stream line by line
  {
    std::istringstream in(line); // make a stream for the line itsel
    int counter2 = 0;
    double element;
    while (in >> element) {
      one_matrix_from_file.at(counter).at(counter2) = element;
      counter2++;
    }
    counter2 = 0;
  }

  // run reduction algorithm for every example and compare to
  for (int i = 0; i < num_examples; i++) {
    Tandem alg(num_orbitals, num_particles, wave_function.at(i));
    vector<double> two_matrix = alg.run();
    Test test(alg);
    vector<vector<double>> one_matrix = test.one_rdm();
    for (size_t j = 0; j < num_orbitals; j++) {
      for (size_t k = 0; k < num_orbitals; k++) {
        BOOST_REQUIRE(abs(one_matrix.at(j).at(k) - one_matrix_from_file.at(j).at(k)) <
		      0.0001);
      }
    }
  }
}
