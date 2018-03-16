/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Sep 30, 2016

  Standalone cpp Poisson solver

*/

//****************** hdf5write.h **************************

#include <string>

#include "H5Cpp.h"
using namespace std;

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

int ReadHDF5File3(string, string, int, int, int, double*); // Rank 3 files
