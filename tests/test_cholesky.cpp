// This program reads .redg density matrices from file and cholesky
// decomoses them generating new .cho file with the same name.
#include <assert.h>
#include <fstream>
#include <iostream> // std::cout
#include <math.h>   // sqrt
#include <random>
#include <string>
#include <vector>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp> // not sure if needed
#include "../cholesky.h"
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE CholeskyDecompTest
#include <boost/test/included/unit_test.hpp>  // include this to get main()

using namespace std;
void print_matrix(vector<double> m) {
  cout << "\n";
  for (size_t i = 0; i < m.size(); i++) {
    cout << m.at(i) << " ";
  }
  cout << "\n";
}

BOOST_AUTO_TEST_CASE( cholesky_decomp_test )
{
  // Ground state to N=3, M=4
  vector<double> two_matrix = {1, 0, 0, 0, 0, 0,
			         1, 0, 0, 0, 0,
			            0, 0, 0, 0,
			               1, 0, 0,
			                  0, 0,
			                     0};
  vector<double> unit_matrix = {1/2.0, 0, 0, 0, 0, 0,
			          1/2.0, 0, 0, 0, 0,
				     1/2.0, 0, 0, 0,
			                1/2.0, 0, 0,
			                   1/2.0, 0,
			                      1/2.0};
  vector<double> res = log_cholesky_decomp(two_matrix, unit_matrix, 4);
  size_t counter = 0;
  vector<double> expected = {-2.5e-08,
			     0, -2.5e-08,
			     0, 0, -8.40562, 
			     0, 0, 0, -2.5e-08,
			     0, 0, 0, 0, -8.40562,
			     0, 0, 0, 0, 0, -8.40562};
  for (size_t i = 0; i < 6; i++) {
    for (size_t j = 0; j <= i; j++) {
      if (abs(exp(expected.at(counter)) - exp(res.at(counter))) > 0.001) {
	cout << "FALIED!\n";
	BOOST_REQUIRE(false);
      }
      counter++;
    }
  }
}
