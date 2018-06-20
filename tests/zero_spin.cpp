#include "../tandem.h"
#include <boost/math/special_functions/binomial.hpp>
#include <vector>
#include "../tandem_function_interface.h"

using namespace std;

int main() {
  cout << "zero_spin distribution test... ";
  // Basic info
  const int num_orbitals = 20;
  const int num_particles = 4;
  const int num_examples = 1;
  const string distribution = "zero_spin";
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  vector<vector<double>> res = tandem(num_particles, num_examples,
				      num_orbitals, distribution);
    cout << " PASSED!" << endl;
  return 0;
}
