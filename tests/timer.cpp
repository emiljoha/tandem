#include "../tandem.h"
#include "../cholesky.h"
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


double abs(vector<double> vec) {
  double res = 0;
  for (size_t i = 0; i < vec.size(); i++) {
    res += vec.at(i) * vec.at(i);
  }
  return sqrt(res);
}
void print_usage(const size_t num_orbitals) {
  cout << "Usage: ./tandem -p num_particles -e num_examples -f file_name -d distriution"
       << endl;
  cout << "Number of orbitals used: " << num_orbitals << endl;
  cout << "Number of orbitals must be known at compile time" << endl;
  cout << "To change set num_orbitals in tandem_cli.cpp:27 and type \"make "
          "cli\" "
       << endl;
  cout << "distribution can be one of \"uniform\",  \"one_over_x-p\" and, \"exponential-T\"" << endl;
  cout << "Where p is the power 1/x^p and T is the temperature exp(-x/T)" << endl;
  return;
}

int main(int argc, char **argv) {
  clock_t tStart = clock();
  cout << "Time to execute tandem... " << endl; 
  size_t num_particles = 0;
  const size_t num_orbitals = 10; // Must be known at compile time
  size_t num_examples = 0;
  string file_name;
  int c;
  std::function<double(int)> distribution_function;
  // Parse commandline options
  opterr = 0;
  if (argc == 1) {
    print_usage(num_orbitals);
    return 1; 
  }
  bool p = false, e = false, f = false, d = false; 
  while ((c = getopt(argc, argv, "p:e:f:d:c")) != -1)
    switch (c) { 
    case 'p':
      num_particles = atoi(optarg);
      p = true; 
      break;
    case 'e':
      num_examples = atoi(optarg);
      e = true;
      break;
    case 'f':
      file_name = optarg;
      f = true;
      break;
    case 'd':
      if (strcoll(optarg, "uniform") == 0){
	distribution_function = [](int n)->double {return 1;}; 
      }
      else if (((string)optarg).substr(0, 10) == "one_over_x"){
	string s(optarg);
	int pos = s.find("-");
	double p = atof(s.substr(pos+1).c_str());
	distribution_function = [p](int n)->double { return 1 / pow(n, p); };
      }
      else if (((string)optarg).substr(0, 11) == "exponential"){
	string s(optarg);
	int pos = s.find("-");
	double T = atof(s.substr(pos+1).c_str());
	distribution_function = [T](int n)->double { return exp(-n / T); };
      }
      else{
	cout << "Unknown distriution " << optarg << endl;
	cout << "Avaliable: \"uniform\" and \"one_over_x\"" << endl;
	cout << "aborting" << endl;
	return -1; 
      }
      d = true;
      break; 
    case '?':
      if (optopt == 'p' or optopt == 'e' or optopt == 'f' or optopt == 'd')
        fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint(optopt))
        fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
      print_usage(num_orbitals); 
      return 1; 
    default:
      abort();
    }
  if (not(p and e and f and d)) {
    cout << "Missing arguments. " << endl; 
    print_usage(num_orbitals);
    return 1;
  }
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  vector<double> wave_function(basis_length);
  ofstream file;
  file.open(file_name);

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
  double element = num_particles * (num_particles - 1) / two_matrix_basis;
  size_t pos = 0;
  for (size_t i = 0; i < two_matrix_basis; i++) {
      previous.at(pos) = element;
      pos = pos + (two_matrix_basis - i);
  };

  // Define a random number generator and initialize it with a reproducible
  // seed.
  // (The seed is unsigned, otherwise the wrong overload may be selected
  // when using mt19937 as the base_generator_type.)
  boost::minstd_rand generator(42);

  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::uniform_real<double> uni_dist(-1,1);
  boost::variate_generator<boost::minstd_rand&, boost::uniform_real<double> > uni(generator, uni_dist);
  // You can now retrieve random numbers from that distribution by means
  // of a STL Generator interface, i.e. calling the generator as a zero-
  // argument function

  for (size_t i = 0; i < num_examples; i++) {
    for (size_t j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = distribution_function(j + 1) * uni();
    }
    double total_sum = abs(wave_function);
    for (size_t j = 0; j < wave_function.size(); j++) {
      wave_function.at(j) = wave_function.at(j) / total_sum;
    }
    assert(abs(abs(wave_function) - 1) < 0.001);
    Tandem alg(num_orbitals, num_particles, wave_function);
    vector<double> two_matrix = alg.run();
    
    vector<double> res;
    string extension = boost::filesystem::extension(file_name);
    if (extension == ".cho"){
      res = log_cholesky_decomp(two_matrix, previous, num_orbitals);
    }
    else if (extension == ".red") {
      res = two_matrix;
    }
    else {
      cout << "Unknown file type " << file_name.substr(pos) << endl;  
      return -1;
    }
    // label
    file << 1 << " "; 
    for (size_t j = 0; j < res.size() - 1; j++) {
      file << res.at(j) << " ";
    }
    file << endl;
    if (i % 10 == 0 or i == num_examples - 1){
      cout << i + 1 << "/" << num_examples << "\r";
      std::cout.flush();
      }
    // Update previous before next loop.
    // previous = two_matrix;
  }
  cout << endl; 
  file.close();
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}
