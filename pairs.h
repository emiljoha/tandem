#ifndef __PAIRS_H_INCLUDED__
#define __PAIRS_H_INCLUDED__

#include <math.h>
#include <vector>
#include <boost/math/special_functions/binomial.hpp> //  boost::math::binomial_coefficient<int>(10, 2);

using namespace std;

bool* get_N_basis_first_state(const int num_orbitals, const int num_particles){
  bool* state = new bool[num_orbitals];
  for (int i = 0; i < num_orbitals; ++i) {
    if (i < num_particles) {
      state[i] = true;
    } else {
      state[i] = false;
    }
  }
  sort(state, state + num_orbitals);
  return state;
}
  
vector<bool> is_spin_zero(const int num_orbitals, const int num_particles,
			  const vector<int> one_basis_spins){
  assert(num_orbitals == one_basis_spins.size());
  // N in the naming stands for N particle state where N is the number of particles.
  size_t N_basis_length = boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  vector<bool> N_basis_spin_is_zero(N_basis_length, false);
  bool* N_basis_state = get_N_basis_first_state(num_orbitals, num_particles);
  size_t num_basis = 0;
  do {
    double total_spin = 0;
    for (int i = 0; i < num_orbitals; i++){
      if (N_basis_state[num_orbitals - 1 - i] == 1){
	total_spin += one_basis_spins.at(i);
      }
    }
    N_basis_spin_is_zero.at(num_basis) = (total_spin == 0);
    total_spin = 0;
    num_basis++;
  } while (next_permutation(N_basis_state, N_basis_state + num_orbitals));
  delete [] N_basis_state;
  return N_basis_spin_is_zero;
}

vector<double> basis_energies(const int num_orbitals, const int num_particles,
			    const vector<double> one_basis_energies){
  assert(num_orbitals == one_basis_energies.size());
  // N in the naming stands for N particle state where N is the number of particles.
  size_t N_basis_length = boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  vector<double> N_basis_energies(N_basis_length, false);
  bool* N_basis_state = get_N_basis_first_state(num_orbitals, num_particles);
  size_t num_basis = 0;
  do {
    double total_energy = 0;
    for (int i = 0; i < num_orbitals; i++){
      if (N_basis_state[num_orbitals - 1 - i] == 1){
	total_energy += one_basis_energies.at(i);
      }
    }
    N_basis_energies.at(num_basis) = total_energy;
    total_energy = 0;
    num_basis++;
  } while (next_permutation(N_basis_state, N_basis_state + num_orbitals));
  delete [] N_basis_state;
  return N_basis_energies;
}
#endif
