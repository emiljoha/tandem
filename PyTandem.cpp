#include "tandem_function_interface.h"
#include "conversions.h"
#include <vector>

using namespace std;
namespace py = boost::python;

BOOST_PYTHON_MODULE(PyTandem) // 
{
  using namespace boost::python;
  // register the vec-to-python converter
  to_python_converter<
    vector<double>,
    std_vector_to_python_list>();
  // register the mat-to-python converter
  to_python_converter<
    vector<vector<double> >,
    std_mat_to_python_list>();
  // register python to vec
  std_vector_from_python_list();
  // Register python to mat
  std_mat_from_python_list();
  def("tandem", tandem,
      return_value_policy<return_by_value>());
}
