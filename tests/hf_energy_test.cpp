// Test that the generated reduced density matrix from the HF D matrix gives the correct energy.
#include "../tandem.h"
#include "../hartree_fock.h"
#include "test.h"
#include <boost/math/special_functions/binomial.hpp>
#include <boost/random/linear_congruential.hpp> // not sure if needed
#include <fstream>
#include <iostream> // std::cout
#include <math.h>   // sqrt
#include <string>
#include <vector>
#include <stdexcept>
#define BOOST_TEST_DYN_LINK        // this is optional
#define BOOST_TEST_MODULE HartreeFockTest
#include <boost/test/included/unit_test.hpp>  // include this to get main()

using namespace std;

vector<vector<double>> loadtxt(string file_name, size_t rows, size_t columns){
  vector<vector<double>> result(rows, vector<double>(columns));
  ifstream file(file_name, ios_base::in);
  if (!file) {
    throw invalid_argument(file_name + " did not open sucessfully ");
  }
  int counter = 0;
  for (std::string line;
       std::getline(file, line);
       counter++) { // read stream line by line
    std::istringstream in(line); // make a stream for the line itself
    int counter2 = 0;
    double element;
    while (in >> element) {
      result.at(counter).at(counter2) = element;
      counter2++;
    }
    counter2 = 0;
  }
  return result;
}

vector<double> vec_loadtxt(string file_name, size_t length) {
  vector<double> result(length);
  ifstream file(file_name, ios_base::in);
  if (!file) {
    throw invalid_argument(file_name + " did not open sucessfully ");
  }
  int counter = 0;
  for (std::string line;
       std::getline(file, line);
       counter++) { // read stream line by line
    if (counter >= 2) {
      throw invalid_argument(file_name + " have multiple lines ");
    }
    std::istringstream in(line); // make a stream for the line itself
    int counter2 = 0;
    double element;
    while (in >> element) {
      result.at(counter2) = element;
      counter2++;
    }
    counter2 = 0;
  }
  return result;
}


double abs(vector<double> vec) {
  double res = 0;
  for (int i = 0; i < vec.size(); i++) {
    res += vec.at(i) * vec.at(i);
  }
  return sqrt(res);
}

vector<vector<double>> transpose(vector<vector<double>> A) {
  vector<vector<double>> result(A.at(0).size(), vector<double>(A.size()));
  for (size_t i = 0; i < A.size(); i++) {
    for (size_t j = 0; j < A.at(0).size(); j++) {
      result.at(j).at(i) = A.at(i).at(j);
    }
  }
  return result;
}

BOOST_AUTO_TEST_CASE( HF_energy_test ) {
  // Basic info
  const int num_orbitals = 20;      // Dependent on hamiltonians!
  const int num_particles = 4;      // Dependent on hamiltonians!
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  // Load HF solution matrix
  vector<vector<vector<double>>> D = {transpose(loadtxt("data/ca_0_D.txt", num_orbitals, num_particles)),
				      transpose(loadtxt("data/ca_0.5_D.txt", num_orbitals, num_particles)),
				      transpose(loadtxt("data/ca_1_D.txt", num_orbitals, num_particles)),
				      transpose(loadtxt("data/ca_1.5_D.txt", num_orbitals, num_particles)),
				      transpose(loadtxt("data/ca_2_D.txt", num_orbitals, num_particles))};

  vector<vector<double>> py_wf = {vec_loadtxt("data/ca_0_wf_py_D.txt", basis_length),
				  vec_loadtxt("data/ca_0.5_wf_py_D.txt", basis_length),
				  vec_loadtxt("data/ca_1_wf_py_D.txt", basis_length),
				  vec_loadtxt("data/ca_1.5_wf_py_D.txt", basis_length),
				  vec_loadtxt("data/ca_2_wf_py_D.txt", basis_length)};

  vector<double> energy_facit = {-34.495999999433018, -35.631316599194399, -36.947231295443920,
				 -38.319448097450618, -39.193570845747900};
  // vector<double> interaction = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2};
  vector<double> interaction = {0, 0.5, 1, 1.5, 2};
  vector<vector<double>> hamiltonian_data = loadtxt("data/H_ca_20_edges.dat", 2, pow(num_orbitals, 4));
  vector<vector<double>> hamiltonian(interaction.size(), vector<double>(pow(num_orbitals, 4)));
  // Magic numbers
  double particle_scale_factor = 1 / (double(num_particles) - 1.0);
  double vrenorm = - pow(42.0 / (40.0 + num_particles), 0.3) / 4.0;
  for (size_t j = 0; j < interaction.size(); j++){
    for (size_t i = 0; i < hamiltonian_data.at(0).size(); i++){
      hamiltonian.at(j).at(i) = -(particle_scale_factor * hamiltonian_data.at(0).at(i) +
				  vrenorm * interaction.at(j) * hamiltonian_data.at(1).at(i));
    }
  }
  // Initialize vectros
  // vector<double> wave_function(basis_length, 0);
  // wave_function.at(0) = 1;
  // cout << "\nHF wave function norm: " << total_sum << endl;
  for (size_t j = 2; j < 3; j++){
    vector<double> wave_function = hf_wf_from_D(D.at(j), num_orbitals, num_particles);
    Tandem alg(num_orbitals, num_particles, wave_function);
    vector<double> two_matrix = alg.run();
    Test test(alg);
    test.expand();
    test.flatten();
    //  = test.two_rdm;
    // vector<vector<double>> one_matrix = test.one_rdm();
    ofstream two_matrix_file("two_matrix_" + to_string(j) + ".txt");
    two_matrix_file.flags (std::ios::scientific);
    two_matrix_file.precision (std::numeric_limits<double>::digits10 + 1);
    // for (size_t i = 0; i < two_matrix.size(); i++){
    //   two_matrix_file << two_matrix.at(i) << " ";
    // }
    for (size_t i = 0; i < num_orbitals; i++){
      for (size_t j2 = 0; j2 < num_orbitals; j2++){
    	for (size_t k = 0; k < num_orbitals; k++){
    	  for (size_t l = 0; l < num_orbitals; l++){
    	    two_matrix_file << test.two_rdm[i][j2][k][l] << " ";
    	  }
    	}
      }
    }
    two_matrix_file.close();
    double total_sum = abs(wave_function);
    Tandem alg_py(num_orbitals, num_particles, py_wf.at(j));
    alg_py.run();
    Test test_py(alg_py);
    // double total_sum = abs(wave_function);
    BOOST_REQUIRE(abs(test_py.energy(hamiltonian.at(j)) - test.energy(hamiltonian.at(j))) < 0.01);
  }
}
