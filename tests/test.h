#ifndef __TEST_H_INCLUDED__
#define __TEST_H_INCLUDED__

#include "../tandem.h"
#include <vector>
#include "boost/multi_array.hpp"

class Test {
public:
  Test(const Tandem &alg); 
  // Functions
  double trace();
  void flatten();
  double energy(vector<double> hamiltonian);
  void expand();
  vector<vector<double> > one_rdm();
  // Variables
  const int num_orbitals;
  const int num_particles;
  boost::multi_array<double, 4> two_rdm;
  bool expanded; 
  bool two_rdm_all_zero;
  bool flattned; 
  bool one_matrix_generated; 
  vector<double> flattned_two_rdm;
  vector<vector<double> > one_matrix;
};

#endif
