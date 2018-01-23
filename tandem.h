#ifndef __TANDEM_H_INCLUDED__
#define __TANDEM_H_INCLUDED__

#include "boost/multi_array.hpp"
#include <vector>

typedef boost::multi_array<double, 4> array_type;
// const auto& comb = boost::math::binomial_coefficient<double>;
using namespace std;

void elementwise_xor(bool* x, bool* y, bool* res, int size);
void elementwise_and(bool* x, bool* y, bool* res, int size);
int count(bool* x, int size);

class Tandem {
public:
  Tandem(const int num_orbitals, const int num_particles, const vector<double> wave_function);
  // Get the 2RDM corresponding to the N-wave_function coefficients
  vector<double> run();
  // Test class need acess to internals
  friend class Test;
private:
  const int num_orbitals;
  const int num_particles;
  bool two_rdm_all_zero;
  const vector<double> wave_function;
  boost::multi_array<double, 4> two_rdm;
  // Private Member Functions
  bool sign(int index[4], bool* state);
  void order_and_save(int index[4], bool* state, int num_state1,
		      int num_state2);
};

#endif

