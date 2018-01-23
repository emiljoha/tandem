// generate tabulated values later used by random_test_cases to test
// consistency of algorithm after changes.

#include "../tandem.h"
#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix.hpp>
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
  cout << "Warning: This will create a new expected behaviur of the algorithm"
       << endl;
  cout << "This should only be run when there is an expected change in the "
          "algorithm"
       << endl;
  cout << "Do you want to continiue? (y/n): ";
  string ans;
  cin >> ans;
  if (ans != "y") {
    cout << endl << "Aborting" << endl;
    return 0;
  }
  cout << "Continuing... " << endl;

  const int num_orbitals = 10;
  const int num_particles = 4;
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  const int num_examples = 10;

  // Set up files to save to
  const string file_name_wf = "random_test_cases.wf";
  const string file_name_dmat = "random_test_cases.red";
  ofstream file_wf(file_name_wf);
  ofstream file_dmat(file_name_dmat);
  // Assert that the files opend succesfully.
  assert(file_wf and file_dmat);

  // Random number generator for the wave functions
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);

  vector<double> wave_function(basis_length, 0);
  // Generate and save num_examples wavefunctions and 2-matrices to files
  for (int i = 0; i < num_examples; i++) {
    for (int j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = distribution(generator);
    }

    double total_sum = abs(wave_function);
    for (int j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = wave_function.at(j) / total_sum;
    }

    Tandem alg(num_orbitals, num_particles, wave_function);
    vector<double> two_matrix = alg.run();
    for (vector<double>::iterator it = two_matrix.begin();
         it != two_matrix.end(); ++it) {
      file_dmat << *it << " ";
    }
    file_dmat << endl;
    for (vector<double>::iterator it = wave_function.begin();
         it != wave_function.end(); ++it) {
      file_wf << *it << " ";
    }
    file_wf << endl;
  }
  file_dmat.close();
  file_wf.close();
  return 0;
}
