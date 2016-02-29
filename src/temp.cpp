/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Dec 3, 2014

  Standalone cpp Poisson solver

*/

//****************** multigrid.cpp **************************

#include "multigrid.h"


void MultiGrid::Trace(Array3D* elec, Array3D** E, double* point, ofstream& file)
{
  // This traces an electron down to the bottom, saving path info if requested
  // Diffusion has now been added. This version recalculates mu at each point.
  // And iterates 1000 steps after reaching the bottom.
  int i, j, k, nsteps = 0, nstepsmax = 10000, bottomsteps = 1000;
  bool ReachedBottom = false;
  double mu, E2, Emag, ve, vth, tau, Tscatt;
  double theta, phiangle, zmin, zbottom;
  double x, y, z, xbot, ybot, zbot;
  zmin = (elec->z[Channelkmax] + elec->z[Channelkmin]) / 2.0;
  zbottom = elec->z[Channelkmin - 1] + elec->dz / 2.0;
  x = point[0]; y = point[1]; z = point[2];
  double*  E_interp = new double[3];
  
  while (nsteps < nstepsmax)    
    {
      nsteps += 1;
      E2 = 0.0;
      for (i=0; i<3; i++)
	{
	  E_interp[i] = E[i]->DataInterpolate3D(point[0],point[1],point[2]);
	  E2 += E_interp[i] * E_interp[i];
	}
      Emag = max(0.1, sqrt(E2));
      mu = mu_Si(Emag * MICRON_PER_CM, CCDTemperature); // Mobility
      ve = mu * MICRON_PER_CM * MICRON_PER_CM; // Drift Velocity Factor (Velocity / E)
      vth = sqrt(3.0 * KBOLTZMANN * CCDTemperature / (2.0 * ME))  * MICRON_PER_M * DiffMultiplier; // Thermal Velocity
      vth = vth / sqrt((double)NumDiffSteps);
      tau  = ME / QE * mu * METER_PER_CM * METER_PER_CM; // scattering time

      phiangle = 2.0 * pi * drand48();
      theta = acos(-1.0 + 2.0 * drand48());
      Tscatt = -tau * log(1.0 - drand48()) * (double)NumDiffSteps;
      point[0] += (vth * sin(theta) * cos(phiangle) + E_interp[0] * ve) * Tscatt;
      point[1] += (vth * sin(theta) * sin(phiangle) + E_interp[1] * ve) * Tscatt;
      point[2] += (vth * cos(theta) + E_interp[2] * ve) * Tscatt;      
      if (point[2] < zmin && !ReachedBottom)
	{
	  ReachedBottom = true;
	  nstepsmax = nsteps + bottomsteps; // After reaching bottom, iterate bottomsteps more steps.
	  //printf("Reached bottom. zmin = %.6f, z = %.6f\n",zmin, point[2]);
	  xbot = 0.0; ybot = 0.0; zbot = 0.0;
	}
      if (ReachedBottom)
	{
	  xbot += point[0]; ybot += point[1]; zbot += point[2];
	  point[2] = max(zbottom, point[2]);
	}
      if (LogPixelPaths == 1)
	{
	  file  << setw(15) << x << setw(15) << y << setw(15) << z << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
	}
    }
  xbot /= (double) bottomsteps; ybot /= (double) bottomsteps; zbot /= (double) bottomsteps;
  // Final location is average of last 1000 steps.
  i = elec->XIndex(xbot);
  j = elec->YIndex(ybot);
  k = elec->ZIndex(zbot);
  //printf("Final Location, z = %.6f, k = %d\n",zbot, k);	      
  //printf("i = %d, j = %d, k = %d\n",i,j,k);
  //fflush(stdout);
  // Note that the requirement k > 0 will delete electrons that reach the last grid cell.
  if (i > 0 && i < elec->nx-1 && j > 0 && j < elec->ny-1 && k > 0 && k < elec->nz-1)
    {
      elec->data[i + j * elec->nx + k * elec->nx * elec->ny] += 1.0;// Add one electron to this grid cell
      //printf("This electron is getting kept\n");
    }
  delete[] E_interp;
  return;
}


void MultiGrid::FindEdge(Array3D* elec, Array3D** E, double* point, double theta, ofstream& file)
{
  sinth = sin(theta);
  costh = cos(theta);
  z0 = point[2];
  x = point[0]; y = point[1];
  pixx = (int)(point[0] / PixelSize);
  pixy = (int)(point[1] / PixelSize);
  lastpixx = pixx; lastpixy = pixy;
  deltar = 1.0;
  tolerance = 0.0001;
  while (deltar > tolerance)
    {
      x += deltar * costh;
      y += deltar * sinth;
      point[0] = x;
      point[1] = y;
      point[2] = z0;
      Trace(elec,E,point,file);      
      newpixx = (int)(point[0] / PixelSize);
      newpixy = (int)(point[1] / PixelSize);
      if (newpixx == pixx && newpixy == pixy)
	{
	  sign = 1.0;
	}
      else
	{
	  sign = -1.0;
	}
      if (newpixx == lastpixx && newpixy == lastpixy)
	{
	  deltar = sign * fabs(deltar);
	}
      else
	{
	  deltar = sign * fabs(deltar) / 2.0;
	}
      lastpixx = newpixx; lastpixy = newpixy;
    }
  point[0] = x;
  point[1] = y;
  point[2] = z0;
  return;
}
