#ifndef __TANDEM_FUNC_INT_INCLUDED__
#define __TANDEM_FUNC_INT_INCLUDED__
#include "tandem_function_interface.h"
#include "tandem.h"
#include "cholesky.h"
#include "pairs.h"
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
#include <assert.h>

double abs(vector<double> vec) {
  double res = 0;
  for (size_t i = 0; i < vec.size(); i++) {
    res += vec.at(i) * vec.at(i);
  }
  return sqrt(res);
}

vector<double> tandem_on_wf(vector<double> wave_function, int num_particles,
			  int num_orbitals){
  // There is a chance that the generated 2-matrix have eigenvalues
  // that are zero. This does not happen very often but when it
  // happens the cholesky decomposition will fail. To remedy this, on
  // those occation we make a convex combination of mostly the actual
  // 2-matrix but also the unit matrix but with trace N(N-1).
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
  Tandem alg(num_orbitals, num_particles, wave_function);
  vector<double> two_matrix = alg.run();
  return log_cholesky_decomp(two_matrix, previous, num_orbitals);
}

vector<vector<double>> tandem_on_wf(vector<vector<double>> wave_functions, int num_particles,
				    int num_orbitals){
  size_t two_matrix_basis = boost::math::binomial_coefficient<double>(num_orbitals, 2);
  size_t two_matrix_length = boost::math::binomial_coefficient<double>(two_matrix_basis, 2)
    + two_matrix_basis;
  vector<vector<double>> results(wave_functions.size(), vector<double>(two_matrix_length));
  for (size_t i = 0; i < wave_functions.size(); i++){
    results.at(i) = tandem_on_wf(wave_functions.at(i), num_particles, num_orbitals);
  }
  return results;
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
  else if (distribution.substr(0, 9) == "zero_spin"){
    // Shame on me, assuming 20 orbitals and Ca40...
    std::cout << "WARNING: zero_spin, distribution is currently hardcoded for a specific Ca40 basis with 20 spin-orbitals." << endl;
    assert(num_orbitals == 20);  // minimum sanity check
    vector<double> energies = {-8.6240, -5.6793, -4.1370, -1.3829};
    vector<int> spins = {7, 3, 1, 5};
    double T = 10.0;
    vector<vector<double>> info = pairs(num_orbitals, num_particles, spins, energies);
    distribution_function = [T, info](int n)->double {
      return info.at(n).at(0) * exp(-(info.at(n).at(1) - info.at(0).at(1)) / T);};
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
  size_t two_matrix_basis = boost::math::binomial_coefficient<double>(num_orbitals, 2);
  size_t two_matrix_length = boost::math::binomial_coefficient<double>(two_matrix_basis, 2);
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
    res.at(i) = tandem_on_wf(wave_function, num_particles, num_orbitals);
    // if (i % 100 == 0){
    //   printf("\r%i/%i", (int)i, num_examples);
    //   fflush(stdout);
    // }
  }
  // printf("\r%i/%i\n", num_examples, num_examples);
  return res;
}

void run_and_save(string file_name, int num_particles, int num_examples,
		  int num_orbitals, string distribution)
{
  vector<vector<double> > res = tandem(num_particles, num_examples,
				       num_orbitals, distribution);
  ofstream file;
  file.open(file_name, std::ios_base::app);
  for (size_t i = 0; i < res.size(); i++) {
    file << 1 << " ";
    for (size_t j = 0; j < res.at(i).size(); j++){
      file << res.at(i).at(j) << " ";
    }
    file << endl;
  }
  file.close();
}

// vector<vector<double> > tandem(int num_particles, int num_examples,
// 			       int num_orbitals, string distribution);
// void run_and_save(string file_name, int num_particles, int num_examples,
// 		  int num_orbitals, string distribution);
#endif
