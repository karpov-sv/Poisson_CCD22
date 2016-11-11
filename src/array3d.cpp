/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Sep 30, 2016

  Standalone cpp Poisson solver

*/

//****************** array3d.cpp **************************

#include "array3d.h"

Array3D::Array3D(double Xmin, double Xmax, int Nx, double Ymin, double Ymax, int Ny, double Zmin, double Zmax, int Nz, double NZExp, double GateOxide)
{
  int i, j, k;
  nx = Nx; ny = Ny; nz = Nz; 
  xmin = Xmin; xmax = Xmax; ymin = Ymin; ymax = Ymax; zmin = Zmin; zmax = Zmax;
  nzexp = NZExp;
  dx = (xmax - xmin) / (double) nx;
  dy = (ymax - ymin) / (double) ny;
  dzp = (zmax - zmin) / (double) nz;
  volume = dx * dy * dzp;
  x = new double[nx]; y = new double[ny]; z = new double[nz]; zp = new double[nz]; zplus = new double[nz]; zminus = new double[nz]; dzpdz = new double[nz]; zpz = new double[nz]; zmz = new double[nz], zw = new double[nz];
  data = new double[nx * ny * nz];

  Channelkmin =   (int)((ZP(GateOxide * EPSILON_SI / EPSILON_OX) - zmin) / dzp) + 1;
  Gox_si = Z(zmin + (double) Channelkmin * dzp, false); // Top of gate oxide, when translated into Si dielectric constant
  Delta_Z = Gox_si * (1.0 - EPSILON_OX / EPSILON_SI);    
  //printf("Nz = %d, Channelkmin = %d, Gox(input) = %.4f, Gox(eff) = %.4f, Gox_si = %4f, Delta_Z = %4f\n",Nz,Channelkmin,GateOxide,Gox_si*EPSILON_OX / EPSILON_SI,Gox_si, Delta_Z);
  // With the non-linear z-axis, it is useful to calculate a number of coordinates
  for (k=0; k<nz; k++)
    {
      zp[k] = zmin + dzp/2.0 + (double) k * dzp; // This is zprime, a uniformly spaced z
      z[k] = Z(zp[k], true); // This is z coordinate corresponding to k
      // zmz[k] is the coordinate of the cell bottom; zpz[k] is the coordinate of the cell top
      // zw[k] is the height of the cell in z
      if (k == 0)
	{
	  zmz[k] = zmin;    
	  zpz[k] = Z(zp[k] + dzp / 2.0, true);
	}
      else if (k == nz - 1)
	{
	  zmz[k] = Z(zp[k] - dzp / 2.0, true);
	  zpz[k] = zmax;
	}
      else
	{
	  zmz[k] = Z(zp[k] - dzp / 2.0, true);
	  zpz[k] = Z(zp[k] + dzp / 2.0, true);	  
	}
      zw[k] = zpz[k] - zmz[k];

      dzpdz[k] = DZPDz(Z(zp[k], false)); // This is the first derivative of z(zp). Used in E calculation
      //printf("k = %d, Z(zp[k]) = %f, dzpdz = %f\n",k,Z(zp[k],false),dzpdz[k]);
      zplus[k] = ZPlus(Z(zp[k],false)); // This combination is used in the Poisson solution
      zminus[k] = ZMinus(Z(zp[k],false));  // This combination is used in the Poisson solution
      if (zplus[k] < 0.0)
	{
	  // This shouldn't happen, but sometimes does for small arrays.
	  // Apparently because of the very coarse estimate of zp"(z).
	  // This is somehwat of a hack, but prevents zplus from going negative.
	  zplus[k] = zminus[k] / 2.0;
	}
      //printf("k = %d, zp = %f, z = %f, zmz[k] = %f, zpz[k] = %f, zw[k] = %f\n",k, zp[k], z[k], zmz[k], zpz[k], zw[k]);      
    }
  /* This may be useful for future unit testsing.
  double ztest;
  int ktest;
  for (k=0; k<100; k++)
    {
      ztest = 0.01 * (double) k;
      ktest = ZIndex(ztest);
      printf("ztest = %f, ktest = %d, z[ktest] = %f\n",ztest, ktest,Z(zp[ktest],true));
    }
  for (k=0; k<10; k++)
    {
      printf("k = %d, ztrue = %f, zfalse = %f, dzpdztrue = %f, dzpdzfalse = %f\n",k,Z(zp[k],true),Z(zp[k],false),DZPDz(Z(zp[k],true)),DZPDz(Z(zp[k],false)));
    }
    exit(0);*/

  for (j=0; j<ny; j++)
    {
      y[j] = ymin + dy/2.0 + (double) j * dy;
    }
  for (i=0; i<nx; i++)
    {
      x[i] = xmin + dx/2.0 + (double) i * dx;	  
      for (j=0; j<ny; j++)
	{
	  for (k=0; k<nz; k++)
	    {
	      data[i + j*nx + k*nx*ny] = 0.0;
	    }
	}
    }
}

Array3D::~Array3D()
{
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] zp;
  delete[] zpz;
  delete[] zmz;
  delete[] zw;  
  delete[] zplus;
  delete[] zminus;
  delete[] dzpdz;
  delete[] data;
}

double Array3D::PyramidalKernel3D(double deltax, double deltay, double deltaz)
{
  if (deltax>=1.0 || deltay>=1.0 || deltaz>=1.0)
    {
      return 0.0;
    }
  else
    {
      return (1.0-deltax)*(1.0-deltay)*(1.0-deltaz);
    }
}

double Array3D::DataInterpolate3D(double xin, double yin, double zin)
{
  int i, j, k, m, n, ml, nl, p, pl;
  double d, norm, PK=0.0, deltax, deltay, deltaz; 
  i=XIndex(xin);
  j=YIndex(yin);
  k=ZIndex(zin);
  d=0.0; norm = 0.0;
  for (m=i-1; m<i+2; m++)
    {
      ml = max(0,min(nx-1,m));
      deltax=fabs((xin-x[ml])/dx);
      for (n=j-1; n<j+2; n++)
	{
	  nl = max(0,min(ny-1,n));
	  deltay=fabs((yin-y[nl])/dy);
	  for (p=k-1; p<k+2; p++)
	    {
	      pl = max(1,min(nz-2,p));
	      if (zin > z[pl]) deltaz = fabs((zin-z[pl])/(z[pl+1]-z[pl]));
	      else deltaz = fabs((zin-z[pl])/(z[pl-1]-z[pl]));	      
	      if (isnan(deltaz))
		{
		  printf("zin = %f, k = %d, pl = %d, zp[pl] = %f, z[pl] = %f, deltaz = %f, PK = %f\n",zin, k, pl, zp[pl], z[pl], deltaz, PK);
		}
    	      PK = PyramidalKernel3D(deltax,deltay,deltaz);

	      norm += PK;
	      d += PK * data[ml + nl * nx + pl * nx * ny];
	    }
	}
    }
  if (norm < 1.0E-15)
    {
      return 0.0;
    }
  else
    {
      return d/norm;
    }
}

double Array3D::ZP(double z)
{
  return - 100.0 * (nzexp - 1.0) * pow(z / 100.0, (nzexp + 1.0)/nzexp) + nzexp * z;
}

double Array3D::DZPDz(double z)
{
  return - (nzexp - 1.0) * (nzexp + 1.0) / nzexp * pow(z / 100.0, 1.0 / nzexp) + nzexp;
}

double Array3D::D2ZPDz2(double z)
{
  return - (nzexp - 1.0) * (nzexp + 1.0) / (nzexp * nzexp) * pow(z / 100.0, 1.0 / nzexp - 1.0) / 100.0;
}

double Array3D::Z(double zp, bool Gox_shift)

// Inverts ZP(z) using Newton's method
// Applies the shift due to the gate oxide if Gox_shift is true.
{
  int i = 0;
  double error = 1.0, lastroot = zp, newroot;
  while (error>1e-12)
    {
      newroot = lastroot - (ZP(lastroot) - zp) / DZPDz(lastroot);
      error = fabs((newroot - lastroot) / lastroot);
      lastroot=newroot;
      i=i+1;
      if (i > 100)
	{
	  printf("Iterations exceeded in Z(zprime). Quitting\n");
	  return newroot;
	}	  
    }
  // Now we apply the shift due to the gate oxide with a different dielectric constant.
  if (Gox_shift)
    {
      if (newroot >= Gox_si) return newroot - Delta_Z;
      else return newroot * EPSILON_OX / EPSILON_SI;  
    }
  else return newroot;
}

double Array3D::ZPlus(double z)
{
  return dx * dx / (dzp * dzp) * (pow(DZPDz(z), 2.0) + dzp / 2.0 * D2ZPDz2(z)); 
}

double Array3D::ZMinus(double z)
{
  return dx * dx / (dzp * dzp) * (pow(DZPDz(z), 2.0) - dzp / 2.0 * D2ZPDz2(z));   
}

int Array3D::XIndex(double x)
{
  return (int)((x - xmin) / dx);
}

int Array3D::YIndex(double y)
{
  return (int)((y - ymin) / dy);
}

int Array3D::ZIndex(double z)
{
    if (z >= Gox_si * EPSILON_OX / EPSILON_SI)
      {
	z += Delta_Z;
      }
    else
      {
	z *= EPSILON_SI / EPSILON_OX;  	                
      }
  return max(0,min(nz-1,(int)((ZP(z) - zmin) / dzp)));
}

