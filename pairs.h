#ifndef __PAIRS_H_INCLUDED__
#define __PAIRS_H_INCLUDED__

#include <math.h>
#include <vector>
#include <boost/math/special_functions/binomial.hpp> //  boost::math::binomial_coefficient<int>(10, 2);

using namespace std;

// returns matrix with dims (num_basis_states, 2). [n][0] 1 if the spin of the n:th basis state is 0 and 0 otherwise. [n][2] is the energy in the one particle case. One particle energies specified with energies argument. Spins vector only need to have the different spins (times 2 to get int) of the shells not all the component. So for example 3 instead of [-3/2, -1/2, 1/2, 3/2]
vector<vector<double>> pairs(const int num_orbitals, const int num_particles, const vector<int> spins, const vector<double> energies)  // For ca40 in 20 orbitals this is spins = [7, 3, 1, 5] and energies = {-8.6240, -5.6793, -4.1370, -1.3829}.
{
  assert(spins.size() == energies.size());
  int orbitals_from_spin = 0;
  for (int spin : spins){
    orbitals_from_spin += spin + 1;
  }
  assert(orbitals_from_spin == num_orbitals);
  size_t basis_length = boost::math::binomial_coefficient<double>(num_orbitals, num_particles);
  vector<vector<double>> res(basis_length, vector<double>(2));
  vector<int> full_spin_list;
  for (int spin : spins){
    for (int m = -spin; m != (spin + 2); m = m + 2){
      full_spin_list.push_back(m);
    }
  }
  vector<double> energy_list;
  for (int i = 0; i < energies.size(); i++){
    for (int m = -spins.at(i); m != (spins.at(i) + 2); m = m + 2){
      energy_list.push_back(energies.at(i));
    }
  }
  // for (int i = full_spin_list.size() - 1; i >= 0; i = i - 1){
  //  cout << abs(full_spin_list.at(i)) << " ";
  // }
  // cout << endl;
  bool* state = new bool[num_orbitals];
  for (int i = 0; i < num_orbitals; ++i) {
    if (i < num_particles) {
      state[i] = 1;
    } else {
      state[i] = 0;
    }
  }
  sort(state, state + num_orbitals);
  size_t counter = 0;
  size_t spin_zero_counter = 0;
  do {
    double total_spin = 0;
    double total_energy = 0; 
    for (int i = 0; i < num_orbitals; i++){
      if (state[num_orbitals - 1 - i] == 1){
	total_spin += full_spin_list.at(i);
	total_energy += energy_list.at(i);
      }
    }
    // spin_zero_counter++;
    // if (spin_zero_counter < 100){
    //   for (int i = 0; i < num_orbitals; i++) {
    // 	cout << state[i] << " ";
    //   }
    //   cout << counter << " " << total_spin << " " << total_energy << endl;
    // }
    if (total_spin == 0) {
      res.at(counter).at(0) = 1;
    }
    else {
      res.at(counter).at(0) = 0;
    }
    res.at(counter).at(1) = total_energy;
    counter++;
  } while (next_permutation(state, state + num_orbitals));
  return res;
}

#endif
