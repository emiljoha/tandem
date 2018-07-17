#ifndef __CONVERSION_H_INCLUDED__
#define __CONVERSION_H_INCLUDED__

/// This is a modified version of this tutorial
/// https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/
#include <boost/python.hpp>
#include <Python.h>

using namespace std;
namespace py = boost::python;

struct std_vector_to_python_list
{
  static PyObject* convert(vector<double> const& v)
      {
	py::list l;
	for (size_t i = 0; i < v.size(); i++) {
	  l.append(v.at(i)); // Well this sucks
	}
        return py::incref(l.ptr());
      }
};

struct std_int_vector_to_python_list
{
  static PyObject* convert(vector<int> const& v)
      {
	py::list l;
	for (size_t i = 0; i < v.size(); i++) {
	  l.append(v.at(i)); // Well this sucks
	}
        return py::incref(l.ptr());
      }
};

struct std_bool_vector_to_python_list
{
  static PyObject* convert(vector<bool> const& v)
      {
	py::list l;
	for (size_t i = 0; i < v.size(); i++) {
	  l.append(v.at(i)); // Well this sucks
	}
        return py::incref(l.ptr());
      }
};

struct std_mat_to_python_list
{
  static PyObject* convert(vector<vector<double>> const& mat)
      {
	py::list l;
	for (size_t i = 0; i < mat.size(); i++) {
	  py::list temp;
	  for (size_t j = 0; j < mat.at(0).size(); j++){
	    temp.append(mat.at(i).at(j)); // Well this sucks
	  }
	  l.append(temp);
	}
        return py::incref(l.ptr());
      }
};

struct std_vector_from_python_list
{
  std_vector_from_python_list()
  {
    // register the vec-to-python converter
    py::converter::registry::push_back(&convertible,
						  &construct,
						  py::type_id<vector<double> >());
  }
  
  // Determine if obj_ptr can be converted in a QString
  static void* convertible(PyObject* obj_ptr)
    {
      if (!PySequence_Check(obj_ptr)) return 0;
      return obj_ptr;
    }

  // Convert obj_ptr into a vector<double>
  static void construct(
    PyObject* obj_ptr,
    py::converter::rvalue_from_python_stage1_data* data)
    {
      // Extract the double data from the python sequence
      vector<double> value(PySequence_Length(obj_ptr));
      for (size_t i = 0; i < value.size(); i++){
	value.at(i) = PyFloat_AsDouble(PySequence_GetItem(obj_ptr, i));
      }
      // Check if all went well
      if (PyErr_Occurred()){
	throw invalid_argument("Cannot convert element to double");
      }
      // Grab pointer to memory into which to construct the new vector
      void* storage = ((py::converter::rvalue_from_python_storage<vector<double> >*)
		       data)->storage.bytes;
 
      // in-place construct the new QString using the character data
      // extraced from the python object
      new (storage) vector<double>(value);
 
      // Stash the memory chunk pointer for later use by boost.python
      data->convertible = storage;
    }
};

struct std_int_vector_from_python_list
{
  std_int_vector_from_python_list()
  {
    // register the vec-to-python converter
    py::converter::registry::push_back(&convertible,
						  &construct,
						  py::type_id<vector<int> >());
  }
  
  // Determine if obj_ptr can be converted in a QString
  static void* convertible(PyObject* obj_ptr)
    {
      if (!PySequence_Check(obj_ptr)) return 0;
      return obj_ptr;
    }

  // Convert obj_ptr into a vector<int>
  static void construct(
    PyObject* obj_ptr,
    py::converter::rvalue_from_python_stage1_data* data)
    {
      // Extract the double data from the python sequence
      vector<int> value(PySequence_Length(obj_ptr));
      for (size_t i = 0; i < value.size(); i++){
	value.at(i) = PyInt_AsLong(PySequence_GetItem(obj_ptr, i));
      }
      // Check if all went well
      if (PyErr_Occurred()){
	throw invalid_argument("Cannot convert element to int");
      }
      // Grab pointer to memory into which to construct the new vector
      void* storage = ((py::converter::rvalue_from_python_storage<vector<int> >*)
		       data)->storage.bytes;
 
      // in-place construct the new QString using the character data
      // extraced from the python object
      new (storage) vector<int>(value);
 
      // Stash the memory chunk pointer for later use by boost.python
      data->convertible = storage;
    }
};

struct std_bool_vector_from_python_list
{
  std_bool_vector_from_python_list()
  {
    // register the vec-to-python converter
    py::converter::registry::push_back(&convertible,
						  &construct,
						  py::type_id<vector<bool> >());
  }
  
  // Determine if obj_ptr can be converted in a QString
  static void* convertible(PyObject* obj_ptr)
    {
      if (!PySequence_Check(obj_ptr)) return 0;
      return obj_ptr;
    }

  // Convert obj_ptr into a vector<bool>
  static void construct(
    PyObject* obj_ptr,
    py::converter::rvalue_from_python_stage1_data* data)
    {
      // Extract the double data from the python sequence
      vector<bool> value(PySequence_Length(obj_ptr));
      for (size_t i = 0; i < value.size(); i++){
	value.at(i) = PyObject_IsTrue(PySequence_GetItem(obj_ptr, i));
      }
      // Check if all went well
      if (PyErr_Occurred()){
	throw invalid_argument("Cannot convert element to bool");
      }
      // Grab pointer to memory into which to construct the new vector
      void* storage = ((py::converter::rvalue_from_python_storage<vector<bool> >*)
		       data)->storage.bytes;
 
      // in-place construct the new QString using the character data
      // extraced from the python object
      new (storage) vector<bool>(value);
 
      // Stash the memory chunk pointer for later use by boost.python
      data->convertible = storage;
    }
};

struct std_mat_from_python_list
{
  std_mat_from_python_list()
  {
    // register the vec-to-python converter
    py::converter::registry::push_back(&convertible,
				       &construct,
				       py::type_id<vector<vector<double> > >());
  }
  
  // Determine if obj_ptr can be converted in a QString
  static void* convertible(PyObject* obj_ptr)
    {
      if (!PySequence_Check(obj_ptr)) return 0;
      if (!PySequence_Check(PySequence_GetItem(obj_ptr, 0))) return 0;
      return obj_ptr;
    }

  // Convert obj_ptr into a vector<vector<double>>
  static void construct(
    PyObject* obj_ptr,
    py::converter::rvalue_from_python_stage1_data* data)
    {
      // Get dimensions of the matrix
      size_t dim1 = PySequence_Length(obj_ptr);
      size_t dim2 = PySequence_Length(PySequence_GetItem(obj_ptr, 0));
      // Extract the double data from the python sequence of sequences
      vector<vector<double> > value(dim1, vector<double>(dim2));
      for (size_t i = 0; i < dim1; i++){
	for (size_t j = 0; j < dim2; j++){
	  value.at(i).at(j) = PyFloat_AsDouble(PySequence_GetItem(PySequence_GetItem(obj_ptr, i), j));
	}	  
      }
      // Check if all went well
      if (PyErr_Occurred()){
	throw invalid_argument("Cannot convert element to double");

      }
      // Grab pointer to memory into which to construct the new vector
      void* storage = ((py::converter::rvalue_from_python_storage<vector<vector<double> > >*)
		       data)->storage.bytes;
 
      // in-place construct the new QString using the character data
      // extraced from the python object
      new (storage) vector<vector<double> > (value);
      // Stash the memory chunk pointer for later use by boost.python
      data->convertible = storage;
    }
};

#endif
