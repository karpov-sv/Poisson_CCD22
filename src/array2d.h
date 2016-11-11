/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Sep 30, 2016

  Standalone cpp Poisson solver

*/

//****************** array2d.h **************************

#include <stdio.h>       
#include <stdlib.h>      
#include <math.h>        

#include <globals.h>


class Array2D //This packages the 2D data sets
{
 public:
  int nx, ny, *data;
  double xmin, xmax, ymin, ymax, dx, dy, dz, *x, *y;
  Array2D() {};
  Array2D(double, double, int, double, double, int);
  ~Array2D();
  int XIndex(double);
  int YIndex(double);

};
