#include "../tandem.h"
#include "test.h"
#include "boost/multi_array.hpp"
#include <boost/math/special_functions/binomial.hpp> //  boost::math::binomial_coefficient<int>(10, 2);
#include <math.h>
#include <vector>

Test::Test(const Tandem &alg)
  :two_rdm(alg.two_rdm), num_orbitals(alg.num_orbitals),
   num_particles(alg.num_particles), expanded(false),
   flattned(false), one_matrix_generated(false),
   one_matrix(num_orbitals, vector<double>(num_orbitals)),
   flattned_two_rdm(num_orbitals * num_orbitals * num_orbitals *
		     num_orbitals) {};

double Test::energy(vector<double> hamiltonian) {
  assert(hamiltonian.size() == flattned_two_rdm.size());
  if (not flattned) {
    flatten();
  }
  double res = 0;
  for (size_t i = 0; i < hamiltonian.size(); i++) {
    res += hamiltonian.at(i) * flattned_two_rdm.at(i);
  }
  return res;
}

void Test::flatten() {
  if (not expanded) {
    expand();
  }
  int counter = 0;
  for (int i = 0; i < num_orbitals; i++) {
    for (int j = 0; j < num_orbitals; j++) {
      for (int k = 0; k < num_orbitals; k++) {
	for (int l = 0; l < num_orbitals; l++) {
	  flattned_two_rdm.at(counter) = two_rdm[i][j][k][l];
	  counter++;
	}
      }
    }
  }
  flattned = true;
}

double Test::trace() {
  expand();
  // trace condition: tr(two_rdm) = num_particles * (num_particles - 1)
  // sum_ij p_ijij
  double trace = 0;
  for (int i = 0; i < num_orbitals; i++) {
    for (int j = 0; j < num_orbitals; j++) {
      trace += two_rdm[i][j][i][j];
    }
  }
  return trace;
}

vector<vector<double> > Test::one_rdm() {
  if (not expanded) {
    expand();
  }
  for (size_t i = 0; i < num_orbitals; i++) {
    for (size_t j = 0; j < num_orbitals; j++) {
      for (size_t k = 0; k < num_orbitals; k++) {
	one_matrix.at(i).at(j) += two_rdm[i][k][j][k] / (num_particles - 1);
      }
    }
  }
  one_matrix_generated = true;
  return one_matrix;
}

void Test::expand(){
    expanded = true;
    // Fill upp density matrix.
    int counter = 0;
    int basis_length =
        boost::math::binomial_coefficient<double>(num_orbitals, 2);
    vector<vector<int>> list_ij(basis_length, vector<int>(2));
    for (int i = 0; i < num_orbitals; i++) {
      for (int j = i + 1; j < num_orbitals; j++) {
        list_ij.at(counter).at(0) = i;
        list_ij.at(counter).at(1) = j;
        counter++;
      }
    }
    for (int n = 0; n < list_ij.size(); n++) {
      for (int m = n; m < list_ij.size(); m++) {
        int i = list_ij.at(n).at(0);
        int j = list_ij.at(n).at(1);
        int k = list_ij.at(m).at(0);
        int l = list_ij.at(m).at(1);
        two_rdm[i][j][l][k] = -two_rdm[i][j][k][l];
        two_rdm[j][i][k][l] = -two_rdm[i][j][k][l];
        two_rdm[j][i][l][k] = two_rdm[i][j][k][l];
        two_rdm[l][k][i][j] = -two_rdm[i][j][k][l];
        two_rdm[l][k][j][i] = two_rdm[i][j][k][l];
        two_rdm[k][l][i][j] = two_rdm[i][j][k][l];
        two_rdm[k][l][j][i] = -two_rdm[i][j][k][l];
      }
    }
  }
