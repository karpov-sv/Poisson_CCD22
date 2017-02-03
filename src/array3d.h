/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Sep 30, 2016

  Standalone cpp Poisson solver

*/

//****************** array3d.h **************************

#include <stdio.h>       
#include <stdlib.h>      
#include <math.h>        

#include <globals.h>


class Array3D //This packages the 3D data sets
{
 public:
  int nx, ny, nz;
  double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dzp, nzexp, sensorthickness, *x, *y, *z, *zp, *zpz, *zmz, *zw, *zplus, *zminus, *dzpdz, *data;
  double Gox_si, Delta_Z;
  int Channelkmin;
  Array3D() {};
  Array3D(double, double, int, double, double, int, double, double, int, double, double, double);
  ~Array3D();
  double PyramidalKernel3D(double, double, double);
  double DataInterpolate3D(double, double, double);
  int XIndex(double);
  int YIndex(double);
  int ZIndex(double);  
  double ZP(double);
  double DZPDz(double);
  double D2ZPDz2(double);
  double Z(double, bool);
  double ZPlus(double);
  double ZMinus(double);
};
