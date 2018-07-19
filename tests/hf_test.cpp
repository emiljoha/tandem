// hartree fock state construction test
#include "../hartree_fock.h"
#include <vector>
#include <iostream>
#include <math.h>
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE HartreeFockTest
#include <boost/test/included/unit_test.hpp>  // include this to get main()

using namespace std;

BOOST_AUTO_TEST_CASE( hartree_fock_test )
{
  size_t num_particles = 2;
  size_t num_orbitals = 3;
  vector<vector<double>> D(num_particles, vector<double>(num_orbitals, 0));
  D.at(0).at(0) = -1 / sqrt(2); D.at(1).at(0) = 1 / sqrt(2);
  D.at(0).at(1) = 1 / sqrt(4); D.at(1).at(1) = 1 / sqrt(4);
  D.at(0).at(2) = 1 / sqrt(4); D.at(1).at(2) = 1 / sqrt(4);
  vector<double> result = hf_wf_from_D(D, num_orbitals, num_particles);
  vector<double> expected_result = {1 / sqrt(2), 1 / sqrt(2), 0};
  BOOST_REQUIRE(result.size() == expected_result.size());
  for (size_t i = 0; i < result.size(); i++){
    BOOST_REQUIRE(abs(result.at(i) - expected_result.at(i)) < 1e-2);
  }
}
