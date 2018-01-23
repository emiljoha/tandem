#include <vector>
#include <boost/math/special_functions/binomial.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <assert.h>
#include <math.h>       /* log */

using namespace std;

void printmat(gsl_matrix * m) {
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      cout << gsl_matrix_get(m, i, j) << " ";
    }
    cout << endl;
  }
  cout << endl;
}

vector<double> log_cholesky_decomp(vector<double> &A, const vector<double> &P, const int num_orbitals){
  gsl_set_error_handler_off ();
  const size_t basis_length =
      boost::math::binomial_coefficient<double>(num_orbitals, 2);
  const size_t num_elements =
    boost::math::binomial_coefficient<double>(basis_length, 2) + basis_length;
  assert(num_elements == A.size());
  gsl_matrix* m = gsl_matrix_calloc(basis_length, basis_length);
  size_t counter = 0;
  for (size_t i = 0; i < basis_length; i++) {
    for (size_t j = i; j < basis_length; j++) {
      gsl_matrix_set(m, j, i, A.at(counter));
      counter++;
    }
  }
  
  if (gsl_linalg_cholesky_decomp(m) == GSL_EDOM) {
    // m is not positive definite but positive semi-definite.  the
    // vector P is the previus vector that is guaranteed to be
    // positive definite.
    assert(A.size() == P.size());
    double pertubation = 1e-6;
    for (size_t i = 0; i < A.size(); i++){
      A.at(i) = (1 - pertubation) * A.at(i) + pertubation * P.at(i);
    };
    gsl_matrix_free(m);
    return log_cholesky_decomp(A, P, num_orbitals);
  }
  else {
    vector<double> res(A.size(), 0);
    counter = 0;
    for (size_t i = 0; i < basis_length; i++) {
      for (size_t j =0; j <= i; j++) {
	if (i == j) {
	  res.at(counter) = log(gsl_matrix_get(m, i, j));
	}
	else {
	  res.at(counter) = gsl_matrix_get(m, i, j);
	}
	counter++;
      }
    }
    gsl_matrix_free(m);
    return res;
  }
}
  
      
