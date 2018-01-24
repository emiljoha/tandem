#include "tandem.h"
#include "cholesky.h"
#include "conversions.h"
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
#include <stdexcept>

using namespace std;
namespace py = boost::python;

double abs(vector<double> vec) {
  double res = 0;
  for (size_t i = 0; i < vec.size(); i++) {
    res += vec.at(i) * vec.at(i);
  }
  return sqrt(res);
}

vector<vector<double> > tandem(int num_particles, int num_examples,
			       int num_orbitals, string distribution)
{
  std::function<double(int)> distribution_function;
  // Create distribution function from string. 
  if (distribution == "uniform"){
    distribution_function = [](int n)->double {return 1;}; 
  }
  else if (distribution.substr(0, 10) == "one_over_x"){
    int pos = distribution.find("-");
    double p = atof(distribution.substr(pos+1).c_str());
    distribution_function = [p](int n)->double { return 1 / pow(n+1, p); };
  }
  else if (distribution.substr(0, 11) == "exponential"){
    int pos = distribution.find("-");
    double T = atof(distribution.substr(pos+1).c_str());
    distribution_function = [T](int n)->double { return exp(-n / T); };
  }
  else{
    cout << "Unknown distriution " << optarg << endl;
    cout << "Avaliable: \"uniform\", \"one_over_x\", and \"exponential\"" << endl;
    cout << "Remember to add number to the ened of one_over_x and exponential like" << endl;
    cout << "one_over_x-2 to have 1/x^2 decay" << endl;
    cout << "aborting" << endl;
    throw std::invalid_argument("Unknown distribution");
  }
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  vector<double> wave_function(basis_length);
  // There is a chance that the generated 2-matrix have eigenvalues
  // that are zero. This does not happen very often but when it
  // happens the cholesky decomposition will fail. To remedy this, on
  // those occation we make a convex combination of mostly the actual
  // 2-matrix but also the previus two matrix. However the first
  // 2-matrix generated will not have a previoius 2-matrix so we
  // initiate that to the unit matrix but with trace N(N-1).
  size_t two_matrix_basis = boost::math::binomial_coefficient<double>(num_orbitals, 2);
  size_t two_matrix_length = boost::math::binomial_coefficient<double>(two_matrix_basis, 2)
    + two_matrix_basis;
  vector<double> previous(two_matrix_length, 0);
  double element = num_particles * (num_particles - 1) / (double) two_matrix_basis;
  size_t pos = 0;
  for (size_t i = 0; i < two_matrix_basis; i++) {
      previous.at(pos) = element;
      pos = pos + (two_matrix_basis - i);
  };

  // Create instance to hold all the glorios 2-RDMs
  vector<vector<double> > res(num_examples, vector<double>(two_matrix_length)); 

  // Define a random number generator and initialize it with a reproducible
  // seed.
  // (The seed is unsigned, otherwise the wrong overload may be selected
  // when using mt19937 as the base_generator_type.)
  boost::minstd_rand generator(std::time(0));

  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::uniform_real<double> uni_dist(-1,1);
  boost::variate_generator<boost::minstd_rand&, boost::uniform_real<double> > uni(generator, uni_dist);
  // You can now retrieve random numbers from that distribution by means
  // of a STL Generator interface, i.e. calling the generator as a zero-
  // argument function

  for (size_t i = 0; i < num_examples; i++) {
    for (size_t j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = distribution_function(j) * uni();
    }
    double total_sum = abs(wave_function);
    for (size_t j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = wave_function.at(j) / total_sum;
    }
    assert(abs(abs(wave_function) - 1) < 0.001);
    Tandem alg(num_orbitals, num_particles, wave_function);
    vector<double> two_matrix = alg.run();
    res.at(i) = log_cholesky_decomp(two_matrix, previous, num_orbitals);
    if (i % 100 == 0){
      printf("\r%i/%i", (int)i, num_examples);
      fflush(stdout);
    }
  }
  printf("\r%i/%i\n", num_examples, num_examples);
  return res;
}

BOOST_PYTHON_MODULE(PyTandem) // 
{
  using namespace boost::python;
  // register the vec-to-python converter
  to_python_converter<
    vector<double>,
    std_vector_to_python_list>();
  // register the mat-to-python converter
  to_python_converter<
    vector<vector<double> >,
    std_mat_to_python_list>();
  // register python to vec
  std_vector_from_python_list();
  // Register python to mat
  std_mat_from_python_list();
  def("tandem", tandem,
      return_value_policy<return_by_value>());
}
