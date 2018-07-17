#include "../tandem.h"
#include <boost/math/special_functions/binomial.hpp>
#include <vector>
#include "../tandem_function_interface.h"
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE ZeroSpinTest
#include <boost/test/included/unit_test.hpp>  // include this to get main()

using namespace std;

BOOST_AUTO_TEST_CASE( zero_spin_test )
{
  const int num_orbitals = 20;
  const int num_particles = 4;
  const int num_examples = 1;
  const string distribution = "zero_spin";
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  vector<vector<double>> res = tandem(num_particles, num_examples,
				      num_orbitals, distribution);
  BOOST_CHECK(true);  // Very basic test: if it did not crash and burn we are happy!
}
