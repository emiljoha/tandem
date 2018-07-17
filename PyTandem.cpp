#include "tandem_function_interface.h"
#include "hartree_fock.h"
#include "pairs.h"
#include "energy.h"
#include "conversions.h"
#include <vector>
#include <string>

using namespace std;
using namespace boost::python;
BOOST_PYTHON_MODULE(PyTandem) // 
{
  docstring_options local_docstring_options(true, true, false);
  // register the vec-to-python converters
  to_python_converter<
    vector<double>,
    std_vector_to_python_list>();
  to_python_converter<
    vector<int>,
    std_int_vector_to_python_list>();
  to_python_converter<
    vector<bool>,
    std_bool_vector_to_python_list>();
  // register the mat-to-python converter
  to_python_converter<
    vector<vector<double> >,
    std_mat_to_python_list>();
  // register python to vecs
  std_vector_from_python_list();
  std_int_vector_from_python_list();
  std_bool_vector_from_python_list();
  // Register python to mat
  std_mat_from_python_list();
  def("tandem", tandem,
      return_value_policy<return_by_value>(),
      "Generate N-Representable density matrices\n\n"
      "Generates N-Representable density matrices by sampling the N-particle wave function\n"
      "space and reducing the wave functions to 2RDMs\n\n"
      "# Arguments:\n"
      "    num_particles (int)\n"
      "    num_examples (int)\n"
      "    num_orbitals (int)\n"
      "    distribution (string) --- Specify from what distribution as a function of N-particle\n"
      "        basis index the N-particle wave functions are to be drawn from.\n"
      "        On of 'uniform', 'one_over_x-P', 'exponential-T'. P and T are float parameters.\n"
      "        'one_over_x-P' samples wave functions coefficient according to 1/x^P. 'exponential-T'\n"
      "        according to e^(-x/T)\n\n"
      "# Returns:\n"
      "    2RDMs (2D list) --- The log-Cholesky decomposition of the 2D matrix representation\n"
      "        of the full two particle reduced density matrix. See Nakata(2017) page 16\n"
      "        for specification of how the rank 4 tensor is transformed to a rank 2 matrix.\n"
      "        See Pinheiro(1996) for details on the log-Cholesky decomposition.\n"
      );
  def<vector<double> (vector<double>, int, int)>("tandem_on_wf", tandem_on_wf,
      return_value_policy<return_by_value>(),
      "Get 2RDM of a wave function\n\n"
      "# Arguments:\n"
      "    wave_function (list) --- The wave function to be reduced\n"
      "    num_particles (int)\n"
      "    num_orbitals (int)\n\n"
      "# Return:\n"
      "    2RDM (list) --- The log-Cholesky decomposition of the 2D matrix representation\n"
      "        of the full two particle reduced density matrix. See Nakata(2017) page 16\n"
      "        for specification of how the rank 4 tensor is transformed to a rank 2 matrix.\n"
      "        See Pinheiro(1996) for details on the log-Cholesky decomposition.\n"
      );
  def<vector<vector<double>> (vector<vector<double>>, int, int)>("tandem_on_wf", tandem_on_wf,
      return_value_policy<return_by_value>(),
      "Get 2RDM of a list of wave functions\n\n"
      "# Arguments:\n"
      "    wave_functions (list) --- The wave functions to be reduced\n"
      "    num_particles (int)\n"
      "    num_orbitals (int)\n\n"
      "# Return:\n"
      "    2RDM (2D list) --- The log-Cholesky decompositions of the 2D matrix representation\n"
      "        of the full two particle reduced density matrix. See Nakata(2017) page 16\n"
      "        for specification of how the rank 4 tensor is transformed to a rank 2 matrix.\n"
      "        See Pinheiro(1996) for details on the log-Cholesky decomposition.\n"
      );
  def("hf_wf_from_D", hf_wf_from_D,
      return_value_policy<return_by_value>(),
      "Get wave function from Hartree-Fock matrix D\n\n"
      "# Arguments:\n"
      "    D (2D list) --- Hartree-Fock basis change matrix. D[particle_number][orbital_number]\n"
      "    num_orbitals (int)\n"
      "    num_particles (int)\n\n"
      "# Returns:\n"
      "    wave_function (list)\n"
      );
  def("energy", calc_energy,
      return_value_policy<return_by_value>(),
      "Get Energy\n\n"
      "Return energy of wave_function with respect to the hamiltonian\n\n"
      "# Arguments:\n"
      "    wave_function (list)\n"
      "    hamiltonian (list)\n"
      "    num_particles (int)\n"
      "    num_orbitals (int)\n\n"
      "# Return:\n"
      "    energy (float)\n");
  def("spin", calc_energy,
      return_value_policy<return_by_value>(),
      "Get Energy\n\n"
      "Return energy of wave_function with respect to the hamiltonian\n\n"
      "# Arguments:\n"
      "    wave_function (list)\n"
      "    hamiltonian (list)\n"
      "    num_particles (int)\n"
      "    num_orbitals (int)\n\n"
      "# Return:\n"
      "    energy (float)\n");
    def("is_spin_zero", is_spin_zero,
      return_value_policy<return_by_value>(),
      "Get list bool wheather each basis state has spin == 0\n\n"
      "# Arguments:\n"
      "    num_orbitals (int)\n\n"
      "    num_particles (int)\n"
      "    one_basis_spins (list) Spin of each one partical orbital basis.\n"
      "# Return:\n"
      "    List of booleans weather each basis has spin == 0 or not.\n");
    def("basis_energies", basis_energies,
      return_value_policy<return_by_value>(),
      "Get list energies of basis states assuming no interaction\n\n"
      "# Arguments:\n"
      "    num_orbitals (int)\n\n"
      "    num_particles (int)\n"
      "    one_basis_energies (list) Spin of each one partical orbital basis.\n"
      "# Return:\n"
      "    List of energies\n");
}
