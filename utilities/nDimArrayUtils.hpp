///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      functions for dynamic allocation and handling of n-dimensional arrays
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#ifndef NDIMARRAYUTILS_HPP
#define NDIMARRAYUTILS_HPP


#include <vector>


#ifdef __CUDACC__  // source is compiled by nvcc
#warning "generating CUDA device functions for nDimArray.hpp"
#define DEVICE __host__ __device__
#else
#define DEVICE
#endif


namespace rpwa {


  // macros for n-dimensional arrays made of nested vectors
#ifndef protect__  // macro that converts arguments containing commas into single string; useful for more complex types T
#define protect__(x...) x
#endif
#define vector2(T) std::vector<std::vector<T> >
#define vector3(T) std::vector<std::vector<std::vector<T> > >
#define vector4(T) std::vector<std::vector<std::vector<std::vector<T> > > >


  // management function for arrays based on pointer^n variables that
  // can be passed to functions without knowing the array size at
  // compile time
  //
  // elements can be accessed in the canonical way using []^n operators
  // this implementation is, due to the levels of indirection,
  // potentially less performant than pseudo-arrays, which on the other
  // hand make element access more cumbersome
  // 
  // array dimensions are passed as unsigned int arrays, where the order
  // reflects the order of the indices:
  //
  // array[x_0][x_1] ... [x_(n - 1)]
  //        |    |            |
  //     dim[0]  |            |
  //          dim[1] ...   dim[n - 1]
  //
  // array memory has to be freed in reverser order of its allocation

  template<typename D, typename T>
  inline
  void
  delete2DArray(D**&    array,   // two-dimensional array to delete
		const T dim[2])  // extents of two-dimensional array
  {
    for (T i = 0; i < dim[0]; ++i)
      if (array[i])
	delete[] array[i];
    delete[] array;
    array = NULL;
  }

  template<typename D, typename T>
  inline
  void
  allocate2DArray(D**&     array,           // two-dimensional array to create
		  const T  dim[2],          // extents of two-dimensional array
		  const D* defaultVal = 0)  // optional default value
  {
    if (array)
      delete2DArray<D, T>(array, dim);
    array = new D*[dim[0]];
    for (T i = 0; i < dim[0]; ++i) {
      array[i] = new D[dim[1]];
      if (defaultVal)
	for (T j = 0; j < dim[1]; ++j)
	  array[i][j] = *defaultVal;
    }
  }


  template<typename D, typename T>
  inline
  void
  delete3DArray(D***&   array,   // three-dimensional array to delete
		const T dim[3])  // extents of three-dimensional array
  {
    for (T i = 0; i < dim[0]; ++i)
      if (array[i])
	delete2DArray<D, T>(array[i], &dim[1]);      
    delete[] array;
  }

  template<typename D, typename T>
  inline
  void
  allocate3DArray(D***&    array,           // three-dimensional array to create
		  const T  dim[3],          // extents of three-dimensional array
		  const D* defaultVal = 0)  // optional default value
  {
    if (array)
      delete3DArray<D, T>(array, dim);
    array = new D**[dim[0]];
    for (T i = 0; i < dim[0]; ++i)
      allocate2DArray<D, T>(array[i], &dim[1], defaultVal);
  }


  // functions for pseudo-n-dimensional arrays where dimensions are
  // mapped onto a one-dimensional array
  //
  // here index is (like in C++ multi-dimensional arrays) row-major, that
  // is values for last index i_(n - 1) are consecutive in memory
  //
  // array[x_0][x_1] ... [x_(n - 1)]
  //        |    |            |
  //     dim[0]  |            |
  //          dim[1] ...   dim[n - 1]
  //
  // no range checks whatsoever are performed

  template<typename D, typename T>
  inline
  T
  allocatePseudoNdimArray(D*&      array,           // pointer to one-dimensional array
			  const T* dim,             // extents of n-dimensional array
			  const T  nmbDim,          // number of dimensions
			  const D* defaultVal = 0)  // optional default value
  {
    T nmbElements = dim[0];
    for (T i = 1; i < nmbDim; ++i)
      nmbElements *= dim[i];
    array = new D[nmbElements];
    if (defaultVal)
      for (T i = 0; i < nmbElements; ++i)
	array[i] = *defaultVal;
    return nmbElements * sizeof(D);
  }


  template<typename D, typename T>
  inline
  T
  allocatePseudoNdimArray(D*&                   array,           // pointer to one-dimensional array
			  const std::vector<T>& dim,             // extents of n-dimensional array
			  const D*              defaultVal = 0)  // optional default value
  {
    return allocatePseudoNdimArray<D, T>(array, &(*(dim.begin())), dim.size(), defaultVal);
  }


  template<typename T>
  inline DEVICE
  T
  indicesToOffset(const T* indices,  // indices to map to one-dimensional array index
		  const T* dim,      // extents of n-dimensional array
		  const T  nmbDim)   // number of dimensions
  {
    T offset = indices[0];
    for (T i = 1; i < nmbDim; ++i)
      offset = offset * dim[i] + indices[i];
    return offset;
  }

  template<typename T>
  inline
  T
  indicesToOffset(const std::vector<T>& indices,  // indices to map to one-dimensional array
		  const std::vector<T>& dim)      // extents of n-dimensional array
  {
    T offset = indices[0];
    for (T i = 1; i < dim.size(); ++i)
      offset = offset * dim[i] + indices[i];
    return offset;
  }


  template<typename T>
  inline DEVICE
  void
  offsetToIndices(const T  offset,   // one-dimensional array index
		  const T* dim,      // extents of n-dimensional array
		  const T  nmbDim,   // number of dimensions
		  T*       indices)  // indices to map onto
  {
    T index = offset;
    for (T i = nmbDim - 1; i >= 1; --i) {
      indices[i] = index % dim[i];
      index      = index / dim[i];
    }
    indices[0] = index;
  }

  template<typename T>
  inline
  void
  offsetToIndices(const T               offset,   // one-dimensional array index
		  const std::vector<T>& dim,      // extents of n-dimensional array
		  std::vector<T>&       indices)  // indices to map onto
  {
    T index = offset;
    for (T i = dim.size() - 1; i >= 1; --i) {
      indices[i] = index % dim[i];
      index      = index / dim[i];
    }
    indices[0] = index;
  }


}  // namespace rpwa


#endif  // NDIMARRAYUTILS_HPP
