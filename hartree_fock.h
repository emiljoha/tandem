#ifndef __HF_H_INCLUDED__
#define __HF_H_INCLUDED__

using namespace std;

#include <vector>
#include <assert.h>
#include <math.h> // sort
#include <boost/math/special_functions/binomial.hpp> //  boost::math::binomial_coefficient
#include <iostream>
int permutation_sign(int* A, int length){
  int result = 1;
  for (int i = 0; i < length; i++){
    for (int j = i + 1; j < length; j++){
      if (A[i] < A[j]) {
	result = -result;
      }
    }
  }
  return result;
}

double coefficient(const int* index, const vector<vector<double>> &D, int num_particles){
  double coefficient = 1;
  for (int i = 0; i < num_particles; i++){
    coefficient = coefficient * D.at(i).at(index[i]);
  }
  return coefficient;
}


vector<double> hf_wf_from_D(vector<vector<double>> D, size_t num_orbitals, size_t num_particles){
  assert(D.size() == num_particles);
  assert(D.at(0).size() == num_orbitals);
  int basis_length = boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  vector<double> result(basis_length);
  int* state = new int[num_orbitals];
  for (int i = 0; i < num_orbitals; ++i) {
    if (i < num_particles) {
      state[i] = 1;
    } else {
      state[i] = 0;
    }
  }
  sort(state, state + num_orbitals);
  size_t state_count = 0;
  int particle_count = 0;
  double total_coefficient = 0;
  do {
    int* index = new int[num_particles];
    particle_count = 0;
    for (int i = num_orbitals - 1; i >= 0; i = i - 1){
      if (state[i] == 1){
	index[particle_count] = num_orbitals - 1 - i;
	particle_count++;
      }
    }
    sort(index, index + num_particles);
    do{
      total_coefficient += coefficient(index, D, num_particles) * permutation_sign(index, num_particles);
    } while (std::next_permutation(index, index + num_particles));
    result.at(state_count) = total_coefficient;
    state_count++;
    total_coefficient = 0;
  } while (std::next_permutation(state, state + num_orbitals));
  return result;
}

#endif
