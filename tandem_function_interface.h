#ifndef __TANDEM_FUNC_INT_INCLUDED__
#define __TANDEM_FUNC_INT_INCLUDED__
vector<vector<double> > tandem(int num_particles, int num_examples,
			       int num_orbitals, string distribution);
void run_and_save(string file_name, int num_particles, int num_examples,
		  int num_orbitals, string distribution);
#endif
