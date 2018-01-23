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
#include "ra.h"
#include "cholesky.h"

using namespace std;
int main() {
  // Basic info
  const int num_orbitals = 10;      // Dependent on hamiltonians!
  const int num_particles = 4;      // Dependent on hamiltonians!
  const int num_matrices = 100; // Dependent on hamiltonians!
  const string file_name_density_matrices = "../../data/dmatC0.redg";
  const int basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, 2);
  const int num_variables = boost::math::binomial_coefficient<double>(basis_length, 2)
    + basis_length;
  //// LOAD DENSITY MATRICES FROM FILE ////
  // Initialize vectors
  
  vector<vector<float>> density_matrices(num_matrices,
					 vector<float>(num_variables, 0));
  // reading in hamitonians to test system (4 particles 10 orbitals)
  ifstream file_load(file_name_density_matrices, ios_base::in);
  // Check if file opended succesfully
  if (!file_load) {
    cout << "ERROR: " << file_name_density_matrices
         << " did not open succesfully" << endl;
    return -1;
  }
  int counter = 0;
  vector<float> energy(num_matrices), interaction(num_matrices);
  for (string line; getline(file_load, line); counter++){ // read stream line by line
    std::istringstream in(line); // make a stream for the line itsel
    in >> energy.at(counter) >> interaction.at(counter);
    int counter2 = 0;
    float element;
    while (in >> element) {
      density_matrices.at(counter).at(counter2) = element;
      counter2++;
    }
    counter2 = 0;
    // As .red does not include the last element for trace
    // conservation we need to recreate that. we can get the diagonal
    // element recursivly from red with the formula:
    // diagonal_i = diagonal_{i-1} + num_basis - i + 1, diagonal_0 = 0
    
    // The sum of the diagonal elements should be:
    // num_particles(num_particles - 1) / 2. The 1/2 is beacuse of
    // that in the full rank 4 tensor there is another diagonal
    // element when switching i<->j and k<->l.
    float fixed_trace = num_particles * (num_particles - 1) / 2;
    double trace = 0;
    size_t diagonal = 0;
    for (int i = 0; i < basis_length; i++) {
      trace += density_matrices.at(counter).at(diagonal);
      diagonal += basis_length - i; 
    }
    density_matrices.at(counter).back() = fixed_trace - trace;
  }
  file_load.close();
  //// PREPARE INITIAL PREVIOUS MATRIX ////
  // in the case of the first matrix having zero valued eigenvalues.
  vector<float> previous(num_variables, 0);
  float element = num_particles * (num_particles - 1) / float(2 * basis_length);
  size_t pos = 0;
  for (int i = 0; i < basis_length; i++) {
      previous.at(pos) = element;
      pos += basis_length - i;
  };

  //// FILE TO SAVE TO ////
  ofstream file_save;
  string file_name_cho = "../../data/dmatC0.chog"; ///file_name_density_matrices;
  file_save.open(file_name_cho);
  if (!file_save) {
    cout << "ERROR: " << file_name_cho
         << " did not open succesfully" << endl;
    return -1;
  }

  /// CALCULATE AND SAVE CHOLESKY DECOMPOSITIONS ////
  vector<float> res(density_matrices.at(0).size());
  for (int i = 0; i < num_matrices; i++) {
    res = log_cholesky_decomp(density_matrices.at(i), previous, num_orbitals);
    // previous = density_matrices.at(i);
    // label
    file_save << energy.at(i) << " " << interaction.at(i) << " "; 
    for (size_t j = 0; j < res.size(); j++) {
      file_save << res.at(j) << " ";
    }
    file_save << endl;
    if (i % 10 == 0 or i == num_matrices - 1){
      cout << i + 1 << "/" << num_matrices << "\r";
      std::cout.flush();
      }
  }
  file_save.close();
  return 0;
}
