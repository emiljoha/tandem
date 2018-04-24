// Takes tabulated wave functions and corresponding 2-matrices and
// test that the algorithm still generates the same 2-matrices as
// before.

#include "../tandem.h"
#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream> // std::cout
#include <math.h>   // sqrt
#include <random>
#include <string>
#include <vector>
//#include <sstream>
#include <fstream>
#include <string>
using namespace std;
typedef boost::numeric::ublas::matrix<double> boost_matrix;
double abs(vector<double> vec) {
  double res = 0;
  for (int i = 0; i < vec.size(); i++) {
    res += vec.at(i) * vec.at(i);
  }
  return sqrt(res);
}
int main() {
  cout << "random test cases... ";
  const int num_orbitals = 10;
  const int num_particles = 4;
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  const int dmat_length =
      boost::math::binomial_coefficient<double>(
          boost::math::binomial_coefficient<double>(num_orbitals, 2), 2) +
      boost::math::binomial_coefficient<double>(num_orbitals, 2);
  const int num_examples = 10;
  vector<vector<double>> wave_function(num_examples,
                                      vector<double>(basis_length, 0));
  vector<vector<double>> expected_two_matrix(num_examples,
                                            vector<double>(dmat_length, 0));
  // Set up files to save to
  const string file_name_wf = "random_test_cases_wavefunction.dat";
  const string file_name_dmat = "random_test_cases_matrices.dat";
  ifstream file_wf(file_name_wf);
  ifstream file_dmat(file_name_dmat);
  // Assert that the files opend succesfully.
  assert(file_wf and file_dmat);
  // read wave function from disk.
  int counter = 0;
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
  // read expected dmat from disk.
  counter = 0;
  for (std::string line;
       std::getline(file_dmat, line) and counter < num_examples;
       counter++) // read stream line by line
  {
    std::istringstream in(line); // make a stream for the line itsel
    int counter2 = 0;
    double element;
    while (in >> element) {
      expected_two_matrix.at(counter).at(counter2) = element;
      counter2++;
    }
    counter2 = 0;
  }
  // Will not need files any more.
  file_dmat.close();
  file_wf.close();
  // Check that files are reasonable
  assert(expected_two_matrix.size() == wave_function.size());
  // Check that Tandem get the same 2-matrix of all wave functions
  for (int i = 0; i < wave_function.size(); i++) {
    Tandem alg(num_orbitals, num_particles, wave_function.at(i));
    vector<double> two_matrix = alg.run();
    assert(two_matrix.size() == expected_two_matrix.at(i).size());
    for (int j = 0; j < two_matrix.size(); j++) {
      if (abs(two_matrix.at(j) - expected_two_matrix.at(i).at(j)) > 0.00001) {
        cout << "ERROR! " << endl;
        cout << "Expected: " << expected_two_matrix.at(i).at(j) << endl;
        cout << "GOT: " << two_matrix.at(j) << endl;
        return 0;
      }
    }
  }
  cout << "PASSED!" << endl;
  return 0;
}
