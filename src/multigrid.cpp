/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Sep 30, 2016

  Standalone cpp Poisson solver

*/

//****************** multigrid.cpp **************************

#include "multigrid.h"

MultiGrid::MultiGrid(string inname) //Constructor
{
  // This reads in the data from the poisson.cfg
  // file, sets the initial conditions, and solves
  // Poisson's equation using SOR and Multi-Grid methods
  double setup_time, solution_time, efield_time, trace_time;
  time_t time1, time2;
  string StepNum;
  string underscore = "_";
  int m, m_init = 0;
  time1 = time(NULL);

  //Set the random number seed
  srand48( (unsigned int) time(NULL));

  // First we read in the configuration information
  ReadConfigurationFile(inname);
  printf("Finished Reading config file\n");
  // Then, we build the multigrid arrays and set the initial conditions
  phi = new Array3D*[nsteps+1];
  rho = new Array3D*[nsteps+1];
  elec = new Array3D*[nsteps+1];
  hole = new Array3D*[nsteps+1];
  eps = new Array3D*[nsteps+1];
  BCType = new Array2D*[nsteps+1];
  QFe = new Array2D*[nsteps+1];
  QFh = new Array2D*[nsteps+1];
  E = new Array3D*[3];

  BuildArrays(phi, rho, elec, hole, eps, E, BCType, QFe, QFh);
  // Save coordinate grid along each axis.
  SaveGrid();
  printf("Finished Building Arrays. \n");
  Setkmins(rho, phi, eps);
  Channelkmin = rho[0]->ChannelCkmin;
  ChannelStopkmin = rho[0]->ChannelStopCkmin;
  SetInitialVoltages(phi[0], BCType[0]);
  //CountCharges(rho, elec, hole);
  SetFixedCharges(rho[0], BCType[0]); // Place fixed charges
  //CountCharges(rho, elec, hole);
  Set_QFh(QFh, qfh); // Set hole Quasi-Fermi level
  time2 = time(NULL);
  setup_time = difftime(time2, time1);
  printf("Finished Setting ICs. Setup time = %.3f seconds.\n",setup_time);
  fflush(stdout);
  if (BuildQFeLookup == 1)
    {
      printf("Building QFe Lookup table\n");
      Adjust_QFe(QFe, rho[0], elec[0]);        
      VCycle_Inner(phi, rho, elec, hole, eps, BCType, QFe, QFh, w, nsteps, ncycle);
      Adjust_QFe(QFe, rho[0], elec[0]);        
      WriteQFeLookup(outputfiledir, outputfilebase, "QFe");
      StepNum = boost::lexical_cast<std::string>(999);      
      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Elec", elec[0]);
      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Hole", hole[0]);
      BuildQFeLookup = 0;
    }
  else
    {
      ReadQFeLookup(outputfiledir, outputfilebase, "QFe");
    }
  // Now we run NumSteps cycle, adding NumElec electrons each step and re-solving
  // Poisson's equation at each step.
  // If a Continuation is requested, we read in the existing results and start
  // from where we left off.
  if (Continuation == 1)
    {
      StepNum = boost::lexical_cast<std::string>(LastContinuationStep);
      ReadOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Elec", elec[0]);
      ReadOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Hole", hole[0]);        
      m_init = LastContinuationStep + 1;
    }

  for (m=m_init; m<NumSteps; m++)
    {
      time1 = time(NULL);
      // Now we cycle through the VCycles to solve Poisson's equation
      for (int n=0; n<iterations; n++)
	{
          printf("--- Starting VCycle %d / %d of step %d / %d.\n",
            n + 1, iterations, m + 1, NumSteps);
	  Adjust_QFe(QFe, rho[0], elec[0]); // Set QFe	  	  
	  if (ElectronAccumulation == 0) FillRho(rho[0], elec[0]); // Add electrons to rho

	  if (m == m_init)
	    {
	      VCycle_Inner(phi, rho, elec, hole, eps, BCType, QFe, QFh, w, nsteps, ncycle);
	    }
	  else
	    {
	      VCycle_Zero(phi, rho, elec, hole, eps, BCType, QFe, QFh, w, ncycle_loop);
	    }
	  SetFixedCharges(rho[0], BCType[0]); // Clear mobile charges out of rho	  
	  //CountCharges(rho, elec, hole);
	}
      StepNum = boost::lexical_cast<std::string>(m);
      time2 = time(NULL);
      solution_time = difftime(time2, time1);
      printf("Finished solving Poisson's equation. Solution time = %.3f seconds\n",solution_time);
      time1 = time(NULL);
      // Next we calculate the E fields
      Gradient(phi[0], E);
      time2 = time(NULL);
      efield_time = difftime(time2, time1);
      printf("Finished calculating E Fields. E Field time = %.3f seconds\n",efield_time);
      // Now, we write out the potential and charge density results
      if (SaveElec !=0 && m % SaveElec == 0)
	{
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Elec", elec[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Hole", hole[0]);
	}
      
      if (SaveData !=0 && m % SaveData == 0)
	{
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "phi", phi[0]);
	  WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "rho", rho[0]);
	  if(LogEField > 0)
	    {
	      
	      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ex", E[0]);
	      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ey", E[1]);
	      WriteOutputFile(outputfiledir, outputfilebase+underscore+StepNum, "Ez", E[2]);
	    }
	}
      time1 = time(NULL);
      // Now we trace the electrons.
      if (PixelBoundaryTestType == 0)
	{
	  TraceGrid(m);
	}
      if (PixelBoundaryTestType == 1)
	{
	  TraceSpot(m);
	  if (ElectronAccumulation == 1) WriteCollectedCharge(outputfiledir, outputfilebase+underscore+StepNum, "CC");
	}
      if (PixelBoundaryTestType == 3)
	{
	  TraceList(m);
	}
      if (PixelBoundaryTestType == 2 || PixelBoundaryTestType == 4)
	{
	  TraceRegion(m);
	}
      // Calculate pixel areas after tracing electrons, if requested.
      if (PixelAreas >= 0 && (m % PixelAreas) == 0)
	{
	  CalculatePixelAreas(m);
	}
      time2 = time(NULL);
      trace_time = difftime(time2, time1);
      printf("Finished tracing electrons. Trace time = %.3f seconds\n",trace_time);
    }
  return;
}

MultiGrid::~MultiGrid() //Destructor                                                                                            
{
  int n;
  for (n=0; n<nsteps+1; n++)
    {
      delete phi[n];
      delete rho[n];
      delete elec[n];
      delete hole[n];
      delete eps[n];            
      delete BCType[n];
      delete QFe[n];
      delete QFh[n];                  
    }
  for (n=0; n<3; n++)
    {
      delete E[n];
    }
  delete[] E;
  delete[] phi;
  delete[] rho;
  delete[] elec;
  delete[] hole;
  delete[] eps;
  delete[] BCType;
  delete[] QFe;
  delete[] QFh;  
  return;
}

void MultiGrid::SaveGrid() {
    printf("Saving coordinate grids.\n");
    string grid_name = outputfiledir + "/grid_";

    // The phi, rho, Ex, Ey, Ez arrays all use the same grid.
    Array3D *A = phi[0];

    string xgrid_name = grid_name + "x.dat";
    ofstream xgrid_out(xgrid_name.c_str());
    xgrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
    for(int i = 0; i < A->nx; ++i) {
      xgrid_out << setw(8) << i << setw(16) << A->x[i] - 0.5 * A->dx << setw(16) << A->x[i] << setw(16) << A->x[i] + 0.5 * A->dx << endl;
    }
    xgrid_out.close();

    string ygrid_name = grid_name + "y.dat";
    ofstream ygrid_out(ygrid_name.c_str());
    ygrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
    for(int i = 0; i < A->ny; ++i) {
      ygrid_out << setw(8) << i << setw(16) << A->y[i] - 0.5 * A->dy << setw(16) << A->y[i] << setw(16) << A->y[i] + 0.5 * A->dy << endl;
    }
    ygrid_out.close();

    string zgrid_name = grid_name + "z.dat";
    ofstream zgrid_out(zgrid_name.c_str());
    zgrid_out << setw(8) << "index" << setw(16) << "lower edge" << setw(16) << "grid center" << setw(16) << "upper edge" << endl;
    for(int i = 0; i < A->nz; ++i) {
      zgrid_out << setw(8) << i << setw(16) << A->zmz[i] << setw(16) << A->z[i] << setw(16) << A->zpz[i] << endl;
    }
    zgrid_out.close();

}


void MultiGrid::ReadPhotonList(string outputfiledir, string name)
{
  // This reads the PhotonList
  int i = 0, success;
  string line, filename = outputfiledir+'/'+name;
  ifstream pl_in(filename.c_str());
  getline (pl_in,line); //Skip first line, which is labels
  while (!pl_in.eof())
    {
      getline (pl_in,line);
      i += 1;
    }
  NumPhotons = i-1;
  pl_in.clear();
  pl_in.seekg (0, ios::beg); //Reset back to line 1
  getline (pl_in,line);//Skip first line, which is labels
  PhotonListx = new double[NumPhotons];
  PhotonListy = new double[NumPhotons];
  PhotonListdxdz = new double[NumPhotons];
  PhotonListdydz = new double[NumPhotons];
  PhotonListlambda = new double[NumPhotons];  
  double x, y, dxdz, dydz, lambda;
  while (!pl_in.eof())
    {
      getline (pl_in,line);
      istringstream iss(line);
      //printf("Line %d = %s\n",i,line.c_str());
      if (!(iss >> i >> x >> y >> dxdz >> dydz >> lambda))
	{
	  success = 0;
	  break; // error
	}
      if (i < 0 || i > NumPhotons - 1)
	{
	  success = 0;
	  break; // error
	}
      PhotonListx[i] = x;
      PhotonListy[i] = y;
      PhotonListdxdz[i] = dxdz;
      PhotonListdydz[i] = dydz;
      PhotonListlambda[i] = lambda;      
    }
  if (i == NumPhotons - 1) success = 1;
  pl_in.close();
  if (success == 1)
    {
      printf("File %s successfully read. Contains %d photons\n", filename.c_str(), NumPhotons);
      return;
    }
  else
    {
      printf("Problem reading file %s. Quitting\n", filename.c_str());
      exit(0);
    }
}

void MultiGrid::WriteQFeLookup(string outputfiledir, string filenamebase, string name)
{
  // This writes the QFeLookup file
  int i;
  double qf, dqf;
  dqf = (QFemax - QFemin) / ((double)NQFe - 1.0);
  string underscore = "_", slash = "/", filename;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ofstream qf_out(filename.c_str());
  qf_out << setw(16) << "QFe" << setw(8) << "Ne" << endl;
  for (i=0; i<NQFe; i++)
    {
      qf = QFemin + dqf * (double)i;
      qf_out << setw(8) << qf << setw(16) << QFeLookup[i] << endl;
    }
  qf_out.close();
  printf("File %s successfully written\n", filename.c_str());
}

void MultiGrid::ReadQFeLookup(string outputfiledir, string filenamebase, string name)
{
  // This reads the QFeLookup file
  int i = 0, success = 0;
  double qfin, dqf, qfcalc, ne;
  dqf = (QFemax - QFemin) / ((double)NQFe - 1.0);  
  string line, underscore = "_", slash = "/", filename;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ifstream qf_in(filename.c_str());
  if (qf_in.is_open())
    {
      getline (qf_in,line); // First line is labels      
      //printf("Label line = %s\n",line.c_str());
      while (!qf_in.eof())
	{
	  getline (qf_in,line);
	  istringstream iss(line);
	  //printf("Line %d = %s\n",i,line.c_str());
	  if (!(iss >> qfin >> ne)) break; // error
	  qfcalc = QFemin + dqf * (double)i;
	  //printf("Line %d. qfin = %f, qfcalc = %f\n",i,qfin, qfcalc);
	  if (fabs(qfcalc - qfin) > 1.0E-4) break; // error	  
	  QFeLookup[i] = ne;
	  i += 1;
	}
      if (i == NQFe) success = 1;
    }
  qf_in.close();
  if (success == 1)
    {
      printf("File %s successfully read\n", filename.c_str());
      return;
    }
  else
    {
      printf("Problem reading file %s. Quitting\n", filename.c_str());
      exit(0);
    }
}


void MultiGrid::WriteCollectedCharge(string outputfiledir, string filenamebase, string name)
{
  // This writes the CollectedCharge file
  int m, q, PixX, PixY;
  string underscore = "_", slash = "/", filename;
  filename = outputfiledir+slash+filenamebase+underscore+name+".dat";
  ofstream cc_out(filename.c_str());
  cc_out << setw(8) << "PixX" << setw(16) << "PixY" << setw(16) << "electrons" << endl;
  for (m=0; m<NumberofPixelRegions; m++)
    {
      for (q=0; q<NumberofFilledWells[m]; q++)
	{
	  // First figure out the pixel coordinates
	  PixX = (int)floor((FilledPixelCoords[m][q][0] - PixelRegionLowerLeft[m][0]) / PixelSize);
	  PixY = (int)floor((FilledPixelCoords[m][q][1] - PixelRegionLowerLeft[m][1]) / PixelSize);
	  cc_out << setw(8) << PixX << setw(16) << PixY << setw(16) << CollectedCharge[m][q] << endl;
	}
    }
  cc_out.close();
  printf("File %s successfully written\n", filename.c_str());
}


void MultiGrid::ReadConfigurationFile(string inname)
{

  qfe = GetDoubleParam(inname, "qfe", 100.0);
  qfh = GetDoubleParam(inname, "qfh", -100.0);
  // Poisson solver constants
  w = GetDoubleParam(inname, "w", 1.9);			// Successive Over-Relaxation factor
  ncycle = GetIntParam(inname, "ncycle", 100);		// Number of SOR cycles at each resolution
  ncycle_loop = GetIntParam(inname, "ncycle_loop", 100);// Number of SOR cycles in repeated loops
  iterations =  GetIntParam(inname, "iterations", 3);	// Number of VCycles
  NZExp = GetDoubleParam(inname, "NZExp", 10.0);        // Non-linear z axis exponent
  
  // Overall setup
  outputfilebase  = GetStringParam(inname,"outputfilebase", "Test"); //Output filename base
  outputfiledir  = GetStringParam(inname,"outputfiledir", "data"); //Output filename directory
  VerboseLevel = GetIntParam(inname, "VerboseLevel", 1); // 0 - minimal output, 1 - normal, 2 - more verbose.
  SaveData =  GetIntParam(inname, "SaveData", 1);     // 0 - Save only Pts, N save phi,rho,E every Nth step
  SaveElec =  GetIntParam(inname, "SaveElec", 1);     // 0 - Save only Pts, N save Elec every Nth step
  ScaleFactor =  GetIntParam(inname, "ScaleFactor", 1);     // Power of 2 that sets the grid size
  // ScaleFactor = 1 means grid size is 5/6 micron, 128 grids in the z-direction
  nsteps = 5 + (int)(log2(ScaleFactor));
  // nsteps is the number of reduction steps in the Vcle_InnerCycle.
  // This gives 4 grid cells in the z-direction at the coarsest scale
  PixelSize = GetDoubleParam(inname, "PixelSize", 10.0);    // Pixel size in microns
  SensorThickness = GetDoubleParam(inname, "SensorThickness", 100.0);    // Sensor Thickness in microns  
  GridsPerPixel = GetIntParam(inname, "GridsPerPixel", 12); // Grids per pixel at ScaleFactor = 1
  GridsPerPixel = GridsPerPixel * ScaleFactor;
  GridSpacing = PixelSize / (double)GridsPerPixel;

  Nx = GetIntParam(inname, "Nx", 160);                // Number of grids in x at ScaleFactor = 1
  Nx = Nx * ScaleFactor;
  Ny = GetIntParam(inname, "Ny", 160);                // Number of grids in y at ScaleFactor = 1
  Ny = Ny * ScaleFactor;
  Nz = GetIntParam(inname, "Nz", 160);                // Number of grids in z at ScaleFactor = 1
  Nz = Nz * ScaleFactor;
  Nzelec = GetIntParam(inname, "Nzelec", 32);                // Number of grids in z in electron and hole grids at ScaleFactor = 1
  Nzelec = Nzelec * ScaleFactor;
  XBCType = GetIntParam(inname, "XBCType", 1);        // 0 - Free BC, 1 - Periodic BC
  YBCType = GetIntParam(inname, "YBCType", 1);        // 0 - Free BC, 1 - Periodic BC
  SimulationRegionLowerLeft = new double[2];
  SimulationRegionLowerLeft[0] = 0.0;  SimulationRegionLowerLeft[1] = 0.0;
  SimulationRegionLowerLeft = GetDoubleList(inname, "SimulationRegionLowerLeft", 2, SimulationRegionLowerLeft);
  // Voltages and Charges
  Vbb = GetDoubleParam(inname, "Vbb", -50.0);		        // Back bias
  Vparallel_lo = GetDoubleParam(inname, "Vparallel_lo", -8.0);	// Parallel Low Voltage
  Vparallel_hi = GetDoubleParam(inname, "Vparallel_hi", 4.0);	// Parallel High Voltage
  Vserial_lo = GetDoubleParam(inname, "Vserial_lo", -6.0);	// Serial Low Voltage
  Vaverage = (8.0 * Vparallel_lo + 4.0 * Vparallel_hi) / 12.0;
  Vchannelstop = GetDoubleParam(inname, "Vchannelstop", 0.0);
  Vscupper = GetDoubleParam(inname, "Vscupper", 5.0);
  GateOxide = GetDoubleParam(inname, "GateOxide", 0.15);
  FieldOxide = GetDoubleParam(inname, "FieldOxide", 0.40);
  FieldOxideTaper = GetDoubleParam(inname, "FieldOxideTaper", 0.50);  
  BackgroundDoping = GetDoubleParam(inname, "BackgroundDoping", -1.0E12);
  ChannelStopDoping = GetDoubleParam(inname, "ChannelStopDoping", -1.0E12);
  ChannelStopDepth = GetDoubleParam(inname, "ChannelStopDepth", 1.0);
  ChannelStopProfile = GetIntParam(inname, "ChannelStopProfile", 0);
  ChannelStopWidth = GetDoubleParam(inname, "ChannelStopWidth", 1.0);
  ChannelStopPeak = GetDoubleParam(inname, "ChannelStopPeak", 0.0);  
  ChannelDoping = GetDoubleParam(inname, "ChannelDoping", -5.0E11);
  ChannelDepth = GetDoubleParam(inname, "ChannelDepth", 1.0);
  ChannelPeak = GetDoubleParam(inname, "ChannelPeak", 0.0);  
  ChannelProfile = GetIntParam(inname, "ChannelProfile", 0);
  TroughImplant = GetIntParam(inname, "TroughImplant", 0);
  TroughDoping = GetDoubleParam(inname, "TroughDoping", 1.0E12);
  TroughProfile = GetIntParam(inname, "TroughProfile", 0);
  TroughDepth = GetDoubleParam(inname, "TroughDepth", 1.0);
  TroughPeak = GetDoubleParam(inname, "TroughPeak", 0.0);  
  TroughWidth = GetDoubleParam(inname, "TroughWidth", 2.0);
  UndepletedChannelStop = GetIntParam(inname, "UndepletedChannelStop", 0);
  HoleConvergenceVoltage = GetDoubleParam(inname, "HoleConvergenceVoltage", 0.1);

  // Continuations
  Continuation = GetIntParam(inname, "Continuation", 0);
  LastContinuationStep = GetIntParam(inname, "LastContinuationStep", 0);    

  // Temperature and diffusion parameters 
  CCDTemperature = GetDoubleParam(inname, "CCDTemperature", 173.0);
  Ni = NS * pow(CCDTemperature / 300.0, 1.5) * exp( - QE * EG / (2.0 * KBOLTZMANN * CCDTemperature)) / pow(MICRON_PER_CM, 3.0);   // Intrinsic carrier concentration per micron^3
  ktq = .026 * CCDTemperature / 300.0;//KBOLTZMANN * CCDTemperature / QE;
  printf("Ni = %g, kT/q = %f\n",Ni,ktq);
  DiffMultiplier = GetDoubleParam(inname, "DiffMultiplier", 1.0);
  SaturationModel = GetIntParam(inname, "SaturationModel", 0);
  NumDiffSteps = GetIntParam(inname, "NumDiffSteps", 1);
  EquilibrateSteps = GetIntParam(inname, "EquilibrateSteps", 100);
  BottomSteps = GetIntParam(inname, "BottomSteps", 1000);
  ElectronAccumulation = GetIntParam(inname, "ElectronAccumulation", 1);  
  NumVertices = GetIntParam(inname,"NumVertices",2);
  ElectronZ0Area = GetDoubleParam(inname,"ElectronZ0Area",100.0);
  ElectronZ0Fill = GetDoubleParam(inname,"ElectronZ0Fill",100.0);
  
  // Pixel Boundary Tests
  bool PixelBoundaryInsideAnyPixelRegion = false, PixelBoundaryInsideThisPixelRegion;
  int i, j, k, jj;
  string regionnum, fillednum;
  LogEField = GetIntParam(inname, "LogEField", 0);
  LogPixelPaths = GetIntParam(inname, "LogPixelPaths", 0);
  PixelAreas = GetIntParam(inname, "PixelAreas", 0);
  PixelBoundaryTestType = GetIntParam(inname, "PixelBoundaryTestType", 0);
  PixelBoundaryLowerLeft = new double(2);
  PixelBoundaryUpperRight = new double(2);
  PixelBoundaryStepSize = new double(2);
  for (j=0; j<2; j++)
    {
      PixelBoundaryLowerLeft[j] = 0.0;
      PixelBoundaryUpperRight[j] = 100.0;
      PixelBoundaryStepSize[j] = 1.0;
    }
  PixelBoundaryLowerLeft = GetDoubleList(inname, "PixelBoundaryLowerLeft", 2, PixelBoundaryLowerLeft);
  PixelBoundaryUpperRight = GetDoubleList(inname, "PixelBoundaryUpperRight", 2, PixelBoundaryUpperRight);
  PixelBoundaryNx = GetIntParam(inname, "PixelBoundaryNx", 9);
  PixelBoundaryNy = GetIntParam(inname, "PixelBoundaryNy", 9);
  NumElec = GetIntParam(inname, "NumElec", 1000);
  NumSteps = GetIntParam(inname, "NumSteps", 100);      

  if (PixelBoundaryTestType == 0)
    {
      PixelBoundaryStepSize = new double(2);
      for (j=0; j<2; j++)
	{
	  PixelBoundaryStepSize[j] = 1.0;
	}
      PixelBoundaryStepSize = GetDoubleList(inname, "PixelBoundaryStepSize", 2, PixelBoundaryStepSize);
    }
  if (PixelBoundaryTestType == 1)
    {
      Sigmax = GetDoubleParam(inname, "Sigmax", 1.0);
      Sigmay = GetDoubleParam(inname, "Sigmay", 1.0);
      Xoffset = GetDoubleParam(inname, "Xoffset", 0.0);
      Yoffset = GetDoubleParam(inname, "Yoffset", 0.0);
    }
  if (PixelBoundaryTestType == 3)
    {
      PhotonList  = GetStringParam(inname,"PhotonList", "None"); //Photon List file
      if (PhotonList == "None")
	{
	  printf("No PhotonList specified.  Quitting.");
	}
      else
	{
	  ReadPhotonList(outputfiledir, PhotonList);
	}
    }

  // Pixel Regions
  NumberofPixelRegions = GetIntParam(inname, "NumberofPixelRegions", 0);
  PixelRegionLowerLeft = new double*[NumberofPixelRegions];
  PixelRegionUpperRight = new double*[NumberofPixelRegions];
  NumberofFilledWells = new int[NumberofPixelRegions];
  CollectingPhases = GetIntParam(inname, "CollectingPhases", 1);
  FilledPixelCoords = new double**[NumberofPixelRegions];
  CollectedCharge = new int*[NumberofPixelRegions];
  ElectronCount = new double*[NumberofPixelRegions];
  //LastElectronCount = new double*[NumberofPixelRegions];
  PixelQFe = new double*[NumberofPixelRegions];
  //LastPixelQFe = new double*[NumberofPixelRegions];

  BuildQFeLookup = GetIntParam(inname, "BuildQFeLookup", 0);
  NQFe = GetIntParam(inname, "NQFe", 81);
  QFemin = GetDoubleParam(inname, "QFemin", 5.0);
  QFemax = GetDoubleParam(inname, "QFemax", 10.0);
  QFeLookup = new double[NQFe];
  
  for (i=0; i<NumberofPixelRegions; i++)
    {
      PixelRegionLowerLeft[i] = new double[2];
      PixelRegionUpperRight[i] = new double[2];
      for (j=0; j<2; j++)
	{
	  PixelRegionLowerLeft[i][j] = 0.0;
	  PixelRegionUpperRight[i][j] = 100.0;
	}
    }
  for (i=0; i<NumberofPixelRegions; i++)
    {
      int NewNumberofFilledWells;
      PixelBoundaryInsideThisPixelRegion = false;
      regionnum = boost::lexical_cast<std::string>(i);
      PixelRegionLowerLeft[i] = GetDoubleList(inname, "PixelRegionLowerLeft_"+regionnum, 2, PixelRegionLowerLeft[i]);
      PixelRegionUpperRight[i] = GetDoubleList(inname, "PixelRegionUpperRight_"+regionnum, 2, PixelRegionUpperRight[i]);
      NumberofFilledWells[i] = GetIntParam(inname, "NumberofFilledWells_"+regionnum, 0);
      NewNumberofFilledWells = NumberofFilledWells[i]; // Will reduce number of filled wells if some are inside PixelBoundary
      if ((PixelBoundaryLowerLeft[0] >= PixelRegionLowerLeft[i][0]) && (PixelBoundaryLowerLeft[1] >= PixelRegionLowerLeft[i][1]) && (PixelBoundaryUpperRight[0] <= PixelRegionUpperRight[i][0]) && (PixelBoundaryUpperRight[1] <= PixelRegionUpperRight[i][1]))
	{
	  PixelBoundaryInsideAnyPixelRegion = true;
	  PixelBoundaryInsideThisPixelRegion = true;
	  NumberofFilledWells[i] += PixelBoundaryNx * PixelBoundaryNy;
	  // PixelBoundary aded to number of wells filled in .cfg file.
	  // Use first PixelBoundaryNx * PixelBoundaryNy locations for PixelBoundary tests,
	  // then any cells filled in .cfg file which are outside pixel boundary
	  NewNumberofFilledWells = NumberofFilledWells[i]; // Will reduce number of filled wells if some are inside PixelBoundary
	  //printf("NewNumberofFilledWells = %d\n",NewNumberofFilledWells);
	}
      CollectedCharge[i] = new int[NumberofFilledWells[i]];
      FilledPixelCoords[i] = new double*[NumberofFilledWells[i]];
      ElectronCount[i] = new double[NumberofFilledWells[i]];
      //LastElectronCount[i] = new double[NumberofFilledWells[i]];
      PixelQFe[i] = new double[NumberofFilledWells[i]];
      //LastPixelQFe[i] = new double[NumberofFilledWells[i]];
      
      for (j=0; j<NumberofFilledWells[i]; j++)
	{
	  int PixX, PixY;
	  PixX = j % PixelBoundaryNx;
	  PixY = (j - PixX) / PixelBoundaryNx;
	  ElectronCount[i][j] = 0.0;
	  //LastElectronCount[i][j] = 0.0;
	  PixelQFe[i][j] = 40.0;
	  //LastPixelQFe[i][j] = 40.0;
	  FilledPixelCoords[i][j] = new double[2];
	  for (k=0; k<2; k++)
	    {
	      FilledPixelCoords[i][j][k] = 0.0;
	    }
	  if (!PixelBoundaryInsideThisPixelRegion)
	    {
	      fillednum = boost::lexical_cast<std::string>(j);
	      CollectedCharge[i][j] = GetIntParam(inname,"CollectedCharge_"+regionnum+"_"+fillednum,0);
	      FilledPixelCoords[i][j] = GetDoubleList(inname, "FilledPixelCoords_"+regionnum+"_"+fillednum, 2, FilledPixelCoords[i][j]);
	    }
	  else
	    {
	      if (j < PixelBoundaryNx * PixelBoundaryNy)
		{
		  CollectedCharge[i][j] = 0.0;
		  FilledPixelCoords[i][j][0] = PixelBoundaryLowerLeft[0] + ((double)PixX + 0.5) * PixelSize;  
		  FilledPixelCoords[i][j][1] = PixelBoundaryLowerLeft[1] + ((double)PixY + 0.5) * PixelSize;   		    
		}
	      else
		{
		  fillednum = boost::lexical_cast<std::string>(j -  PixelBoundaryNx * PixelBoundaryNy);
		  CollectedCharge[i][j] = GetIntParam(inname,"CollectedCharge_"+regionnum+"_"+fillednum,0);
		  FilledPixelCoords[i][j] = GetDoubleList(inname, "FilledPixelCoords_"+regionnum+"_"+fillednum, 2, FilledPixelCoords[i][j]);
		  // Check if this is inside the PixelRegion
		  if ((FilledPixelCoords[i][j][0] >= PixelBoundaryLowerLeft[0]) && (FilledPixelCoords[i][j][1] >= PixelBoundaryLowerLeft[1]) && (FilledPixelCoords[i][j][0] <= PixelBoundaryUpperRight[0]) && (FilledPixelCoords[i][j][1] <= PixelBoundaryUpperRight[1]))
		    {
		      // This pixel is inside the pixel boundary.  Delete it from the list and move the charge to the appropriate pixel
		      NewNumberofFilledWells -= 1;
		      PixX = (int)floor((FilledPixelCoords[i][j][0] - PixelBoundaryLowerLeft[0]) / PixelSize);
		      PixY = (int)floor((FilledPixelCoords[i][j][1] - PixelBoundaryLowerLeft[1]) / PixelSize);
		      jj = PixX + PixelBoundaryNx * PixY;
		      CollectedCharge[i][jj] = CollectedCharge[i][j];
		      CollectedCharge[i][j] = 0.0;
		      printf("PixX = %d, PixY = %d, jj = %d, NewNumberofFilledWells = %d\n",PixX, PixY, jj, NewNumberofFilledWells);
		    }
		}
	    }
	}
      NumberofFilledWells[i] = NewNumberofFilledWells;
    }
  if(!PixelBoundaryInsideAnyPixelRegion)
    {
      printf("Pixel Boundary not Inside any Pixel Region!!!  Quitting\n");
      exit(0);
    }
  // Fixed Voltage Regions
  NumberofFixedRegions = GetIntParam(inname, "NumberofFixedRegions", 0);
  FixedRegionLowerLeft = new double*[NumberofFixedRegions];
  FixedRegionUpperRight = new double*[NumberofFixedRegions];
  FixedRegionVoltage = new double[NumberofFixedRegions];
  FixedRegionDoping = new double[NumberofFixedRegions];
  FixedRegionChargeDepth = new double[NumberofFixedRegions];
  FixedRegionBCType = new int[NumberofFixedRegions];

  for (i=0; i<NumberofFixedRegions; i++)
    {
      FixedRegionLowerLeft[i] = new double[2];
      FixedRegionUpperRight[i] = new double[2];
      for (j=0; j<2; j++)
	{
	  FixedRegionLowerLeft[i][j] = 0.0;
	  FixedRegionUpperRight[i][j] = 100.0;
	}
    }
  for (i=0; i<NumberofFixedRegions; i++)
    {
      regionnum = boost::lexical_cast<std::string>(i);
      FixedRegionLowerLeft[i] = GetDoubleList(inname, "FixedRegionLowerLeft_"+regionnum, 2, FixedRegionLowerLeft[i]);
      FixedRegionUpperRight[i] = GetDoubleList(inname, "FixedRegionUpperRight_"+regionnum, 2, FixedRegionUpperRight[i]);
      FixedRegionVoltage[i] = GetDoubleParam(inname, "FixedRegionVoltage_"+regionnum,0.0);
      FixedRegionDoping[i] = GetDoubleParam(inname, "FixedRegionDoping_"+regionnum,0.0);
      FixedRegionChargeDepth[i] = GetDoubleParam(inname, "FixedRegionChargeDepth_"+regionnum,0.0);
      FixedRegionBCType[i] = GetIntParam(inname, "FixedRegionBCType_"+regionnum,0);
    }


    // Filter band configuration.
    FilterBand = GetStringParam(inname, "FilterBand", "none");
    if(FilterBand == "u") {
        FilterIndex = 0;
    }
    else if(FilterBand == "g") {
        FilterIndex = 1;
    }
    else if(FilterBand == "r") {
        FilterIndex = 2;
    }
    else if(FilterBand == "i") {
        FilterIndex = 3;
    }
    else if(FilterBand == "z") {
        FilterIndex = 4;
    }
    else if(FilterBand == "y") {
        FilterIndex = 5;
    }
    else {
        printf("No filter response will be used.\n");
        FilterIndex = -1;
    }
    FilterFile = GetStringParam(inname, "FilterFile", "notebooks/depth_pdf.dat");
    if(FilterIndex >= 0) {
        ifstream filter_input(FilterFile.c_str());
        string header;
        getline(filter_input, header); // Skip header line
        for(int i = 0; i < n_filter_cdf; i++) {
            for(int j = 0; j < n_band; j++) {
                assert(filter_input >> filter_cdf[j * n_filter_cdf + i]);
            }
        }
        filter_input.close();
    }
  return;
}

void MultiGrid::BuildArrays(Array3D** phi, Array3D** rho, Array3D** elec, Array3D** hole, Array3D** eps, Array3D** E, Array2D** BCType, Array2D** QFe, Array2D** QFh)
{
  // Builds the multigrid arrays
  int nx, ny, nz, nze, nxx, nyy, nzz, n;
  double xmin, xmax, ymin, ymax, zmin, zmax, zmaxelec, dx, dy, dz;
  nxx = Nx + 1;
  nyy = Ny + 1;
  nzz = Nz + 1;
  dx = PixelSize / (double)GridsPerPixel;
  dy = dx;
  dz = SensorThickness / (double)Nz;
  Xmin = SimulationRegionLowerLeft[0];
  Ymin = SimulationRegionLowerLeft[1];
  Zmin = -dz / 2.0;
  Xmax = dx * (double)nxx + Xmin;
  Ymax = dy * (double)nyy + Ymin;
  Zmax = SensorThickness + dz / 2.0;
  phi[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz,NZExp,SensorThickness);
  rho[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz,NZExp,SensorThickness);
  zmaxelec = rho[0]->Z(rho[0]->zp[Nzelec] + rho[0]->dzp / 2.0);
  elec[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,zmaxelec,Nzelec,NZExp,SensorThickness);
  hole[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,zmaxelec,Nzelec,NZExp,SensorThickness);
  eps[0] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,zmaxelec,Nzelec,NZExp,SensorThickness);
  BCType[0] = new Array2D(Xmin,Xmax,nxx,Ymin,Ymax,nyy);
  QFe[0] = new Array2D(Xmin,Xmax,nxx,Ymin,Ymax,nyy);
  QFh[0] = new Array2D(Xmin,Xmax,nxx,Ymin,Ymax,nyy);  

  for (n=1; n<nsteps+1; n++)
    {
      nx = (phi[0]->nx - 1) / (int)pow(2,n) + 1;
      ny = (phi[0]->ny - 1) / (int)pow(2,n) + 1;
      nz = (phi[0]->nz - 1) / (int)pow(2,n) + 1;
      nze = Nzelec/ (int)pow(2,n);      
      dx = phi[0]->dx * (int)pow(2,n);
      dy = phi[0]->dy * (int)pow(2,n);
      dz = phi[0]->dzp * (int)pow(2,n);
      xmin = phi[0]->xmin + phi[0]->dx / 2.0 - dx / 2.0;
      ymin = phi[0]->ymin + phi[0]->dy / 2.0 - dy / 2.0;
      zmin = phi[0]->zmin + phi[0]->dzp / 2.0 - dz / 2.0;
      xmax = phi[0]->xmax - phi[0]->dx / 2.0 + dx / 2.0;
      ymax = phi[0]->ymax - phi[0]->dy / 2.0 + dy / 2.0;
      zmax = phi[0]->zmax - phi[0]->dzp / 2.0 + dz / 2.0;
      phi[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,NZExp,SensorThickness);
      rho[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,NZExp,SensorThickness);
      elec[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nze,NZExp,SensorThickness);
      hole[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nze,NZExp,SensorThickness);
      eps[n] = new Array3D(xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nze,NZExp,SensorThickness);      
      BCType[n] = new Array2D(xmin,xmax,nx,ymin,ymax,ny);
      QFe[n] = new Array2D(xmin,xmax,nx,ymin,ymax,ny);
      QFh[n] = new Array2D(xmin,xmax,nx,ymin,ymax,ny);            
    }

  for (n=0; n<3; n++)
    {
      E[n] = new Array3D(Xmin,Xmax,nxx,Ymin,Ymax,nyy,Zmin,Zmax,nzz,NZExp,SensorThickness);
    }
  return;
}

void MultiGrid::SetInitialVoltages(Array3D* phi, Array2D* BCType)
{
  int i, j, k, n, index, index2;
  int PixX, PixY;
  double PixXmin, PixYmin;

  // Potential on top
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  index = i + j * phi->nx + (phi->nz - 1) * phi->nx * phi->ny;
	  phi->data[index] = Vbb;
	}
    }
    // Fixed Potentials or free Boundary Conditions on bottom
  for (n=0; n<NumberofFixedRegions; n++)
    {
      for (i=0; i<phi->nx; i++)
	{
	  for (j=0; j<phi->ny; j++)
	    {
	      index = i + j * phi->nx;
	      if (phi->x[i] >= FixedRegionLowerLeft[n][0] && phi->x[i] <= FixedRegionUpperRight[n][0] && phi->y[j] >= FixedRegionLowerLeft[n][1] && phi->y[j] <= FixedRegionUpperRight[n][1])
		{
		  if (FixedRegionBCType[n] == 0) //  Fixed potential region
		    {
		      phi->data[index] = FixedRegionVoltage[n];
		    }
		  else // Free Boundary Condition region
		    {
		      BCType->data[index] = 1.0;
		    }
		}
	    }
	}
    }
    // Potentials in Pixel Region
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<phi->nx; i++)
	{
	  if (phi->x[i] < PixelRegionLowerLeft[n][0] || phi->x[i] > PixelRegionUpperRight[n][0])
	    {
	      continue; // If not in PixelRegion, continue
	    }
	  for (j=0; j<phi->ny; j++)
	    {
	      if (phi->y[j] < PixelRegionLowerLeft[n][1] || phi->y[j] > PixelRegionUpperRight[n][1])
		{
		  continue; // If not in PixelRegion, continue
		}
	      index = i + j * phi->nx;
	      PixX = (int)floor((phi->x[i] - PixelRegionLowerLeft[n][0]) / PixelSize);
	      PixY = (int)floor((phi->y[j] - PixelRegionLowerLeft[n][1]) / PixelSize);
	      PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSize;
	      PixYmin = PixelRegionLowerLeft[n][1] + (double)PixY * PixelSize;
	      // Set the gate voltages
	      for (k=0; k<phi->Vkmin[i]; k++)
		{
		  index2 = index + k * phi->nx * phi->ny;
		  if (CollectingPhases == 2) // Two collecting phases
		    {
		      if (phi->y[j] >= PixYmin + PixelSize/6.0 && phi->y[j] <= PixYmin + 5.0*PixelSize/6.0)
			{
			  // This is the collection region
			  phi->data[index2] = Vparallel_hi;
			}
		      else
			{
			  // This is the barrier gate region
			  phi->data[index2] = Vparallel_lo;
			}
		    }
		  else if (CollectingPhases == 1) // One collecting phase
		    {
		      if (phi->y[j] >= PixYmin + PixelSize/3.0 && phi->y[j] <= PixYmin + 2.0*PixelSize/3.0)
			{
			  // This is the collection region
			  phi->data[index2] = Vparallel_hi;
			}
		      else
			{
			  // This is the barrier gate region
			  phi->data[index2] = Vparallel_lo;
			}
		    }
		}
	    }
	}
    }
  printf("Finished setting Boundary Potentials, \n");
  fflush(stdout);
  return;
}

void MultiGrid::SetFixedCharges(Array3D* rho, Array2D* BCType)
{
  int i, j, k, n, index, index2;
  int FixedRegionkmin, FixedRegionkmax;
  int Ckmin = Channelkmin, Tkmin = Troughkmin;
  Gox_effective = (rho->zmz[rho->ChannelCkmin] - rho->z[rho->ChannelVkmin - 1]);
  double Fox_effective = (rho->zmz[rho->ChannelStopCkmin] - rho->z[rho->ChannelStopVkmin - 1]);
  FixedRegionkmin = Channelkmin;
  Troughkmin = Channelkmin;
  if (ChannelProfile == 0) // Square well
    {
      Ckmin = rho->ZIndex(rho->z[Channelkmin] + ChannelPeak);      
      Channelkmax = rho->ZIndex(rho->z[Channelkmin] + ChannelDepth + ChannelPeak);
    }
  if (ChannelProfile == 1) // Gaussian
    {
      Channelkmax = rho->ZIndex(rho->z[Channelkmin] + 4.0 * ChannelDepth + ChannelPeak);
    }
  if (TroughImplant == 1 && TroughProfile == 0) // Square well
    {
      Tkmin = rho->ZIndex(rho->z[Troughkmin] + TroughPeak);                  
      Troughkmax = rho->ZIndex(rho->z[Troughkmin] + TroughDepth + TroughPeak);
    }
  if (TroughImplant == 1 && TroughProfile == 1) // Gaussian
    {
      Troughkmax = rho->ZIndex(rho->z[Troughkmin] + 4.0 * TroughDepth + TroughPeak);
    }
  int CSkmax=0, PixX;
  double PixXmin, FRChargeDepth, CChargeDepth, TChargeDepth;
  double ChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / pow(MICRON_PER_CM, 3);
  // ChargeFactor converts doping in cm^-3 into the appropriate units
  double ChannelCharge, TroughCharge, FixedRegionCharge;
  ChannelCharge =  ChannelDoping * MICRON_PER_CM  * ChargeFactor;
  TroughCharge =  TroughDoping * MICRON_PER_CM  * ChargeFactor;
  CChargeDepth = rho->zpz[Channelkmax] - rho->zmz[Ckmin];
  TChargeDepth = rho->zpz[Troughkmax] - rho->zmz[Tkmin];
  //if(VerboseLevel > 1) printf("CChargeDepth = %.4f\n",CChargeDepth);
  double ChannelZmin, ChannelZmax, ChannelStopZmin, ChannelStopZmax, ChannelTotal, ChannelStopTotal;
  double TroughZmin, TroughZmax, TroughTotal;
  ChannelZmin = rho->zmz[Ckmin];
  ChannelZmax = rho->zpz[Channelkmax];
  ChannelTotal = erf((ChannelZmax - ChannelZmin - ChannelPeak) / (sqrt(2.0) * ChannelDepth)) + erf(ChannelPeak / (sqrt(2.0) * ChannelDepth));
  TroughZmin = rho->zmz[Troughkmin];
  TroughZmax = rho->zpz[Troughkmax];
  TroughTotal = erf((TroughZmax - TroughZmin - TroughPeak) / (sqrt(2.0) * TroughDepth)) + erf(TroughPeak / (sqrt(2.0) * TroughDepth));
  // Set the background charge:
  for (i=0; i<rho->nx; i++)
    {
      for (j=0; j<rho->ny; j++)
		{
		  for (k=rho->Ckmin[i]; k<rho->nz-1; k++)
			{
			  index = i + j * rho->nx + k * rho->nx * rho->ny;
			  rho->data[index] = BackgroundDoping * ChargeFactor;
			}
		}
    }

  // Fixed Potentials or free Boundary Conditions on bottom
  for (n=0; n<NumberofFixedRegions; n++)
    {
      FixedRegionkmax = rho->ZIndex(rho->Z(rho->zp[FixedRegionkmin]) + FixedRegionChargeDepth[n]);
      FRChargeDepth = rho->zpz[FixedRegionkmax] - rho->zmz[FixedRegionkmin];
      FixedRegionCharge =  FixedRegionDoping[i] * MICRON_PER_CM / FRChargeDepth * ChargeFactor;
      for (i=0; i<rho->nx; i++)
		{
		  for (j=0; j<rho->ny; j++)
			{
			  index = i + j * rho->nx;
			  if (rho->x[i] >= FixedRegionLowerLeft[n][0] && rho->x[i] <= FixedRegionUpperRight[n][0] && rho->y[j] >= FixedRegionLowerLeft[n][1] && rho->y[j] <= FixedRegionUpperRight[n][1])
				{
				  if (FixedRegionBCType[n] == 0) //  Fixed potential region
					{
					  for (k=FixedRegionkmin; k<FixedRegionkmax+1; k++)
						{
						  index2 = index + k * rho->nx * rho->ny;
						  rho->data[index2] = FixedRegionCharge;
						}
					}
				  else // Free Boundary Condition region
					{
					  BCType->data[index] = 1.0;
					}
				}
			}
		}
    }

  // Charges in Pixel Region
  ChannelStopkmax = -100;
  for (n=0; n<NumberofPixelRegions; n++)
    {
      for (i=0; i<rho->nx; i++)
		{
		  PixX = (int)floor((rho->x[i] - PixelRegionLowerLeft[n][0]) / PixelSize);
		  PixXmin = PixelRegionLowerLeft[n][0] + (double)PixX * PixelSize;
		  if (rho->x[i] < PixelRegionLowerLeft[n][0] || rho->x[i] > PixelRegionUpperRight[n][0])
			{
			  continue; // If not in PixelRegion, continue
			}
		  if (ChannelStopProfile == 0) // Square well
		    {
		      CSkmax = rho->ZIndex(rho->z[rho->Ckmin[i]] + ChannelStopDepth);
		    }
		  if (ChannelStopProfile == 1) // Gaussian
		    {
		      CSkmax = rho->ZIndex(rho->z[rho->Ckmin[i]] + 1.5 * ChannelStopDepth + ChannelStopPeak);
		    }
		  if (CSkmax > ChannelStopkmax) ChannelStopkmax = CSkmax;
		  ChannelStopZmin = rho->zmz[rho->Ckmin[i]];
		  ChannelStopZmax = rho->zpz[CSkmax];
		  ChannelStopTotal = erf((ChannelStopZmax - ChannelStopZmin - ChannelStopPeak) / (sqrt(2.0) * ChannelStopDepth)) + erf(ChannelStopPeak / (sqrt(2.0) * ChannelStopDepth));
		  CSChargeDepth = ChannelStopZmax - ChannelStopZmin;
		  for (j=0; j<rho->ny; j++)
			{
			  if (rho->y[j] < PixelRegionLowerLeft[n][1] || rho->y[j] > PixelRegionUpperRight[n][1])
				{
				  continue; // If not in PixelRegion, continue
				}
			  index = i + j * rho->nx;
			  // Now set the charges
			  if (rho->x[i] <= PixXmin + ChannelStopWidth/2.0 + FieldOxideTaper|| rho->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0 - FieldOxideTaper)
				{
				  // This is the Channel Stop Region
				  ChannelStopCharge = ChannelStopDoping * MICRON_PER_CM  * ChargeFactor * rho->TR[i];
				  // Reduce charge by the taper ratio to account for outdiffusion.
				  for (k=0; k<rho->Ckmin[i]; k++)
				    {
				      index2 = index + k * rho->nx * rho->ny;
				      rho->data[index2] = 0.0;

				      //if (j==112 && (i>=104 && i<=120))printf("i = %d, j = %d, k = %d, rho = %f\n",i,j,k,rho->data[index2]);
				    }
				  for (k=rho->Ckmin[i]; k<CSkmax+1; k++)
					{
					  index2 = index + k * rho->nx * rho->ny;
					  if (ChannelStopProfile == 0) // Square Well
						{
						  rho->data[index2] = ChannelStopCharge / CSChargeDepth;
						}
					  if (ChannelStopProfile == 1) // Gaussian
						{
						  rho->data[index2] = ChannelStopCharge / ChannelStopTotal / rho->zw[k] * (erf((rho->zpz[k] - ChannelStopZmin - ChannelStopPeak) / (sqrt(2.0) * ChannelStopDepth)) - erf((rho->zmz[k] - ChannelStopZmin - ChannelStopPeak) / (sqrt(2.0) * ChannelStopDepth)));
						}

					  //if (j==112 && (i>=104 && i<=120))printf("i = %d, j = %d, k = %d, rho = %f\n",i,j,k,rho->data[index2]);
					}
				}
			  else
				{
				  // This is the channel region
				  for (k=rho->Ckmin[i]; k<Channelkmax+1; k++)
					{
					  index2 = index + k * rho->nx * rho->ny;
					  if (ChannelProfile == 0) // Square Well
						{
						  rho->data[index2] = ChannelCharge / CChargeDepth;
						}
					  if (ChannelProfile == 1) // Gaussian
						{
						  rho->data[index2] = ChannelCharge / ChannelTotal / rho->zw[k] * (erf((rho->zpz[k] - ChannelZmin - ChannelPeak) / (sqrt(2.0) * ChannelDepth)) - erf((rho->zmz[k] - ChannelZmin - ChannelPeak) / (sqrt(2.0) * ChannelDepth)));						  
						}
					  if (TroughImplant == 1 && rho->x[i] >= PixXmin + PixelSize/2.0 - TroughWidth / 2.0 && rho->x[i] <= PixXmin + PixelSize/2.0 + TroughWidth / 2.0)
					    {
					      // This is the trough region
					      if (TroughProfile == 0) // Square Well
						{
						  rho->data[index2] = TroughCharge / TChargeDepth;
						}
					      if (TroughProfile == 1) // Gaussian
						{
						  rho->data[index2] = TroughCharge / TroughTotal / rho->zw[k] * (erf((rho->zpz[k] - TroughZmin - TroughPeak) / (sqrt(2.0) * TroughDepth)) - erf((rho->zmz[k] - TroughZmin - TroughPeak) / (sqrt(2.0) * TroughDepth)));
						}
					    }
					  //if (j==112 && (i>=104 && i<=120))printf("i = %d, j = %d, k = %d, rho = %f\n",i,j,k,rho->data[index2]);					  

					  
					}
				}
			}
		}
    }
  if(VerboseLevel > 1)
    {
      printf("In SetFixedCharges, Channelkmin = %d, Channelkmax = %d, ChannelZWidth = %d grid cells, ChannelCharge = %.4f \n",Channelkmin, Channelkmax, (Channelkmax - Channelkmin + 1), ChannelCharge);
      printf("In SetFixedCharges, ChannelStopkmin = %d, ChannelStopkmax = %d, ChannelStopZWidth = %d grid cells, ChannelStopCharge = %.4f \n",ChannelStopkmin, ChannelStopkmax, (ChannelStopkmax - ChannelStopkmin + 1), ChannelStopCharge);
      printf("Effective Gate Oxide thickness = %.3f microns, Effective Field Oxide thickness = %.3f microns.\n",Gox_effective, Fox_effective);
    }
  return;
}


double MultiGrid::SOR_Inner(Array3D* phi, Array3D* rho, Array3D* elec, Array3D* hole, Array3D* eps, Array2D* BCType, Array2D* QFe, Array2D* QFh, double w)
{
  // This is the main loop that calculates the potentials in the array through Successive Over Relaxation.
  // An inner Newton's method loop is used to solve the non-linear Quasi-Fermi level equations.
  // Assumes fixed potentials on the top, mixture of fixed and free BC on bottom, free or periodic BC on the sides.
  double newphi, oldnewphi, tol, exponent;
  double omw, w6, hsquared;
  double apx, amx, apy, amy, apz, amz; // Dielectric constant averages
  double SORChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI));
  double ElecCharge, HoleCharge, Term1, Term2, CellVolume, AveIterations = 0.0;
  double NumElec, NumHoles, TotalHoles=0.0, TotalElectrons=0.0;
  int nn = 0, mm = 0, iter_counter, iter_limit = 10000, red_black;
  int i, j, k, kstart, im, ip, j0, jm, jp, nxy, ind, ind2, indmx, indpx, indmy, indpy, indmz, indpz;
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  omw = 1.0 - w;
  for (red_black=0; red_black<2; red_black++)
    {
      for (i=0; i<phi->nx; i++)
	{
	  if (XBCType == 0)
	    {
	      if (i == 0) {im = 0;} else {im = 1;} // Free BC
	      if (i == phi->nx-1) {ip = 0;} else {ip = 1;} // Free BC
	    }
	  else
	    {
	      if (i == 0) {im = -phi->nx + 2;} else {im = 1;} // Periodic BC
	      if (i == phi->nx-1) {ip = -phi->nx + 2;} else {ip = 1;} // Periodic BC
	    }

	  for (j=0; j<phi->ny; j++)
	    {
	      j0 = j * phi->nx;
	      if (YBCType == 0)
		{
		  if (j == 0) {jm = 0;} else {jm = phi->nx;} // Free BC
		  if (j == phi->ny-1) {jp = 0;} else {jp = phi->nx;} // Free BC
		}
	      else
		{
		  if (j == 0) {jm = (-phi->ny + 2) * phi->nx;} else {jm = phi->nx;} // Periodic BC
		  if (j == phi->ny-1) {jp = (-phi->ny + 2) * phi->nx;} else {jp = phi->nx;} // Periodic BC
		}
	      ind2 = i + j0;
	      if (BCType->data[ind2] > 0) // Free BC at z = 0
		{
		  // TODO fix this for epsilons
		  w6 = w / 6.0;		  
		  indmx = ind2 - im;
		  indpx = ind2 + ip;
		  indmy = ind2 - jm;
		  indpy = ind2 + jp;
		  indmz = ind2;
		  indpz = ind2 + nxy;
		  newphi = omw * phi->data[ind2] + w6 * (phi->data[indmx]+phi->data[indpx]+phi->data[indmy]+phi->data[indpy]+phi->data[indmz]+phi->data[indpz] + hsquared * rho->data[ind2]);
		  phi->data[ind2] = newphi;
		}
	      for (k=phi->Vkmin[i]; k<elec->nz; k++)
		{
		  if ((i + j + k + red_black) % 2 == 0) continue; // Implements Red-Black alternation		  
		  ind = ind2 + k * nxy;
		  indmx = ind - im;
		  indpx = ind + ip;
		  indmy = ind - jm;
		  indpy = ind + jp;
		  indmz = ind - nxy;
		  indpz = ind + nxy;
		  if (k < eps->nz - 1)
		    {
		      apx = (eps->data[ind] + eps->data[indpx]) / 2.0;
		      amx = (eps->data[ind] + eps->data[indmx]) / 2.0;	      
		      apy = (eps->data[ind] + eps->data[indpy]) / 2.0;
		      amy = (eps->data[ind] + eps->data[indmy]) / 2.0;	      
		      apz = (eps->data[ind] + eps->data[indpz]) / 2.0;
		      amz = (eps->data[ind] + eps->data[indmz]) / 2.0;	      
		    }
		  else
		    {
		      apx = 1.0; amx = 1.0; apy = 1.0; amy = 1.0; apz = 1.0; amz = 1.0;
		    }
		  w6 = w / (apx + amx + apy + amy + apz * phi->zplus[k] + amz * phi->zminus[k]);
		  CellVolume = rho->dx * rho->dy * rho->zw[k];
		  Term1 = omw * phi->data[ind] + w6 * (amx * phi->data[indmx] + apx * phi->data[indpx] + amy * phi->data[indmy] + apy * phi->data[indpy] + amz * phi->zminus[k] * phi->data[indmz] + apz * phi->zplus[k] * phi->data[indpz] + hsquared * rho->data[ind]);
		  if (k>=phi->Ckmin[i] && phi->data[ind]>QFe->data[ind2]-30.0*ktq)
		    {
		      // Inner Newton loop for electrons
		      mm++;
		      iter_counter = 0;
		      newphi = phi->data[ind];
		      oldnewphi = -1000.0;
		      tol = 1.0;
		      while (tol > 1.0E-9 && iter_counter < iter_limit)
			{
			  oldnewphi = newphi;
			  exponent = min(50.0, (newphi - QFe->data[ind2]) / ktq);// Prevents ElecCharge getting too large
			  ElecCharge = - Ni * exp(exponent) * SORChargeFactor;
			  Term2 = hsquared * w6 * ElecCharge / ktq;
			  if (fabs(1.0 - Term2) < 1.0E-18) break; 
			  newphi = newphi - (newphi - (w6 * hsquared * ElecCharge + Term1)) / (1.0 - Term2);
			  tol = fabs(newphi - oldnewphi);		      
			  iter_counter++;
			}
		      nn += iter_counter;
		      NumElec = -ElecCharge / SORChargeFactor * CellVolume;
		      elec->data[ind] = NumElec;
		      TotalElectrons += NumElec;
		      if (isnan(newphi))
			{
			  printf("Nan encountered in SOR_Inner electrons at point i,j,k = (%d,%d,%d)\n",i,j,k);
			  exit(0);
			}
		      if (iter_counter >= iter_limit)
			{
			  printf("Warning electron inner loop iteration exceeded at point i,j,k = (%d,%d,%d)\n",i,j,k);
			}
		    }
		  else if (k>=phi->Ckmin[i] && phi->data[ind]<QFh->data[ind2]+30.0*ktq)
		    {
		      // Inner Newton loop for holes
		      mm++;
		      iter_counter = 0;
		      newphi = phi->data[ind];
		      oldnewphi = -1000.0;
		      tol = 1.0;
		      while (tol > 1.0E-9 && iter_counter < iter_limit)
			{
			  oldnewphi = newphi;
			  exponent = min(50.0,(QFh->data[ind2] - newphi) / ktq);// Prevents HoleCharge getting too large
			  HoleCharge = Ni * exp(exponent) * SORChargeFactor;
			  Term2 = hsquared * w6 * HoleCharge / ktq;
			  if (fabs(1.0 + Term2) < 1.0E-18) break;
			  newphi = newphi - (newphi - (w6 * hsquared * HoleCharge + Term1)) / (1.0 + Term2);	  
			  tol = fabs(newphi - oldnewphi);		      
			  iter_counter++;
			}
		      nn += iter_counter;
		      NumHoles = HoleCharge / SORChargeFactor * CellVolume;
		      hole->data[ind] = NumHoles;
		      TotalHoles += NumHoles;
		      if (isnan(newphi))
			{
			  printf("Nan encountered in SOR_Inner holes at point i,j,k = (%d,%d,%d)\n",i,j,k);
			  exit(0);
			}
		      if (iter_counter >= iter_limit)
			{
			  printf("Warning hole inner loop iteration exceeded at point i,j,k = (%d,%d,%d)\n",i,j,k);
			  exit(0);
			}
		      //if (i==120 && j==112 && phi->nz==161)printf("(i,j,k)=(%d,%d,%d), phi=%f, newphi=%f, HoleCharge=%g, exponent=%f,ktq=%f\n",i,j,k,phi->data[ind],newphi,HoleCharge,exponent,ktq);
		    }
		  else
		    {
		      // No inner Newton loop if no free carriers.
		      newphi = Term1;
		      if (ElectronAccumulation == 1)
			{
			  elec->data[ind] = 0.0;
			}
		      hole->data[ind] = 0.0;
		    }
		  phi->data[ind] = newphi;
		}
	      kstart = max(elec->nz, phi->Vkmin[i]);
	      for (k=kstart; k<phi->nz-1; k++)
		{
		  // Above the elec and hole arrays, no inner loop is needed
		  if ((i + j + k + red_black) % 2 == 0) continue; // Implements Red-Black alternation		  
		  ind = ind2 + k * nxy;
		  indmx = ind - im;
		  indpx = ind + ip;
		  indmy = ind - jm;
		  indpy = ind + jp;
		  indmz = ind - nxy;
		  indpz = ind + nxy;
		  w6 = w / (4.0 + phi->zplus[k] + phi->zminus[k]);
		  newphi = omw * phi->data[ind] + w6 * (phi->data[indmx] + phi->data[indpx] + phi->data[indmy] + phi->data[indpy] + phi->zminus[k] * phi->data[indmz] + phi->zplus[k] * phi->data[indpz] + hsquared * rho->data[ind]);
		  phi->data[ind] = newphi;
		}
	    }
	}
    }
  if (rho->nz == 321) printf("Finished SOR_Inner, nz = %d, TotalElectrons = %.1f, TotalHoles = %.1f\n",rho->nz,TotalElectrons,TotalHoles);
  if (mm > 0) AveIterations = (double)nn / (double)mm;
  return AveIterations;
}


double MultiGrid::Error_Inner(Array3D* phi, Array3D* rho, Array3D* elec, Array3D* hole, Array3D* eps)
{
  double RhoChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / (rho->dx * rho->dy);
  double newphi, rhosum, error = 0.0;
  double hsquared, w6;
  double apx, amx, apy, amy, apz, amz; // Dielectric constant averages  
  int i, j, k, nxy, ind, indmx, indpx, indmy, indpy, indmz, indpz;
  nxy = phi->nx * phi->ny;
  hsquared =  phi->dx * phi->dy;
  for (i=1; i<phi->nx-1; i++)
    { 
      for (j=1; j<phi->ny-1; j++)
	{
	  for (k=phi->Vkmin[i]; k<phi->nz-1; k++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      indmx = ind - 1;
	      indpx = ind + 1;
	      indmy = ind - phi->nx;
	      indpy = ind + phi->nx;
	      indmz = ind - nxy;
	      indpz = ind + nxy;
	      if (k < elec->nz)
		{
		  rhosum = rho->data[ind] + (hole->data[ind] - elec->data[ind]) * RhoChargeFactor / rho->zw[k];
		}
	      else
		{
		  rhosum = rho->data[ind];
		}
	      if (k < eps->nz - 1)
		{
		  apx = (eps->data[ind] + eps->data[indpx]) / 2.0;
		  amx = (eps->data[ind] + eps->data[indmx]) / 2.0;	      
		  apy = (eps->data[ind] + eps->data[indpy]) / 2.0;
		  amy = (eps->data[ind] + eps->data[indmy]) / 2.0;	      
		  apz = (eps->data[ind] + eps->data[indpz]) / 2.0;
		  amz = (eps->data[ind] + eps->data[indmz]) / 2.0;	      
		}
	      else
		{
		  apx = 1.0; amx = 1.0; apy = 1.0; amy = 1.0; apz = 1.0; amz = 1.0;
		}
	      w6 = 1.0 / (apx + amx + apy + amy + apz * phi->zplus[k] + amz * phi->zminus[k]);
	      newphi = w6 * (amx * phi->data[indmx] + apx * phi->data[indpx] + amy * phi->data[indmy] + apy * phi->data[indpy] + amz * phi->zminus[k] * phi->data[indmz] + apz * phi->zplus[k] * phi->data[indpz] + hsquared * rhosum);
	      error = max(error, fabs(phi->data[ind] - newphi));
	    }
	}
    }
  //if(VerboseLevel > 1) printf("Grid Size = %d, error = %.3g\n",phi->nx, error);
  return error;
}

void MultiGrid::Restrict(Array3D* phi, Array3D* newphi, Array3D* rho, Array3D* newrho, Array2D* BCType, Array2D* newBCType)
{
  // Assumes fixed potentials on the top, mixture of fixed and free BC on bottom, free or periodic BC on the sides.
  int i, j, k, im, i0, ip, jm, j0, jp, km, k0, kp, nxy, newnxy;
  int newindex, ind000;
  int ind00p, ind0p0, indp00;
  int ind00m, ind0m0, indm00;
  int ind0mm, ind0mp, ind0pm, ind0pp;
  int indm0m, indm0p, indp0m, indp0p;
  int indmm0, indmp0, indpm0, indpp0;
  int indpmm, indpmp, indppm, indppp;
  int indmmm, indmmp, indmpm, indmpp;
  double phisum, rhosum;
  nxy = phi->nx * phi->ny;
  newnxy = newphi->nx * newphi->ny;

  for (i=0; i<newphi->nx; i++)
    {
      if (XBCType == 0) // Free BC
	{
	  if (i == 0) {i0 = 0; im = 0; ip = 1;}
	  else if (i == newphi->nx - 1) {i0 = phi->nx - 1; im = phi->nx - 2; ip = phi->nx - 1;}
	  else {i0 = 2 * i; im = 2 * i - 1; ip = 2 * i + 1;}
	}
      else // Periodic BC
	{
	  if (i == 0) {i0 = 0; im = phi->nx - 2; ip = 1;}
	  else if (i == newphi->nx - 1) {i0 = phi->nx - 1; im = phi->nx - 2; ip = 1;}
	  else {i0 = 2 * i; im = 2 * i - 1; ip = 2 * i + 1;}
	}
      for (j=0; j<newphi->ny; j++)
	{
	  if (YBCType == 0) // Free BC
	    {
	      if (j == 0) {j0 = 0; jm = 0; jp = 1;}
	      else if (j == newphi->ny - 1) {j0 = phi->ny - 1; jm = phi->ny - 2; jp = phi->ny - 1;}
	      else {j0 = 2 * j; jm = 2 * j - 1; jp = 2 * j + 1;}
	    }
	  else // Periodic BC
	    {
	      if (j == 0) {j0 = 0; jm = phi->ny - 2; jp = 1;}
	      else if (j == newphi->ny - 1) {j0 = phi->ny - 1; jm = phi->ny - 2; jp = 1;}
	      else {j0 = 2 * j; jm = 2 * j - 1; jp = 2 * j + 1;}
	    }
	  newBCType->data[i + j * newphi->nx] = BCType->data[i0 + j0 * phi->nx];
	  ind000 = i0 + j0 * phi->nx;
	  indm00 = im + j0 * phi->nx;
	  indp00 = ip + j0 * phi->nx;
	  ind0m0 = i0 + jm * phi->nx;
	  ind0p0 = i0 + jp * phi->nx;
	  indmm0 = im + jm * phi->nx;
	  indmp0 = im + jp * phi->nx;
	  indpm0 = ip + jm * phi->nx;
	  indpp0 = ip + jp * phi->nx;
	  for (k=0; k<newphi->nz; k++)
	    {
	      if (k == 0) {k0 = 0; km = 0; kp = 0;}
	      else if (k == newphi->nz - 1) {k0 = phi->nz - 1; km = phi->nz - 1; kp = phi->nz - 1;}
	      else {k0 = 2 * k; km = 2 * k - 1; kp = 2 * k + 1;}
	      newindex = i + j*newphi->nx + k*newnxy;
	      ind000 = i0 + j0 * phi->nx + k0 * nxy;
	      phisum = 8.0 * phi->data[ind000];
	      rhosum = 8.0 * rho->data[ind000];
	      indm00 = im + j0 * phi->nx + k0 * nxy;
	      indp00 = ip + j0 * phi->nx + k0 * nxy;
	      ind0m0 = i0 + jm * phi->nx + k0 * nxy;
	      ind0p0 = i0 + jp * phi->nx + k0 * nxy;
	      ind00m = i0 + j0 * phi->nx + km * nxy;
	      ind00p = i0 + j0 * phi->nx + kp * nxy;
	      phisum += 4.0 * (phi->data[indm00] + phi->data[indp00] + phi->data[ind0m0] + phi->data[ind0p0] + phi->data[ind00m] + phi->data[ind00p]);
	      rhosum += 4.0 * (rho->data[indm00] + rho->data[indp00] + rho->data[ind0m0] + rho->data[ind0p0] + rho->data[ind00m] + rho->data[ind00p]);
	      indmm0 = im + jm * phi->nx + k0 * nxy;
	      indmp0 = im + jp * phi->nx + k0 * nxy;
	      indpm0 = ip + jm * phi->nx + k0 * nxy;
	      indpp0 = ip + jp * phi->nx + k0 * nxy;
	      phisum += 2.0 * (phi->data[indmm0] + phi->data[indmp0] + phi->data[indpm0] + phi->data[indpp0]);
	      rhosum += 2.0 * (rho->data[indmm0] + rho->data[indmp0] + rho->data[indpm0] + rho->data[indpp0]);
	      indm0m = im + j0 * phi->nx + km * nxy;
	      indm0p = im + j0 * phi->nx + kp * nxy;
	      indp0m = ip + j0 * phi->nx + km * nxy;
	      indp0p = ip + j0 * phi->nx + kp * nxy;
	      phisum += 2.0 * (phi->data[indm0m] + phi->data[indm0p] + phi->data[indp0m] + phi->data[indp0p]);
	      rhosum += 2.0 * (rho->data[indm0m] + rho->data[indm0p] + rho->data[indp0m] + rho->data[indp0p]);
	      ind0mm = i0 + jm * phi->nx + km * nxy;
	      ind0mp = i0 + jm * phi->nx + kp * nxy;
	      ind0pm = i0 + jp * phi->nx + km * nxy;
	      ind0pp = i0 + jp * phi->nx + kp * nxy;
	      phisum += 2.0 * (phi->data[ind0mm] + phi->data[ind0mp] + phi->data[ind0pm] + phi->data[ind0pp]);
	      rhosum += 2.0 * (rho->data[ind0mm] + rho->data[ind0mp] + rho->data[ind0pm] + rho->data[ind0pp]);
	      indmmm = im + jm * phi->nx + km * nxy;
	      indmmp = im + jm * phi->nx + kp * nxy;
	      indmpm = im + jp * phi->nx + km * nxy;
	      indmpp = im + jp * phi->nx + kp * nxy;
	      indpmm = ip + jm * phi->nx + km * nxy;
	      indpmp = ip + jm * phi->nx + kp * nxy;
	      indppm = ip + jp * phi->nx + km * nxy;
	      indppp = ip + jp * phi->nx + kp * nxy;
	      phisum += phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp];
	      rhosum += rho->data[indmmm] + rho->data[indmmp] + rho->data[indmpm] + rho->data[indmpp] + rho->data[indpmm] + rho->data[indpmp] + rho->data[indppm] + rho->data[indppp];
	      newphi->data[newindex] = phisum / 64.0;
	      newrho->data[newindex] = rhosum / 64.0;
	    }
	}
    }
  return;
  }



void MultiGrid::Prolongate(Array3D* phi, Array3D* newphi, Array2D* newBCType)
{
  // Assumes fixed potentials on the top, mixture of fixed and free BC on bottom, free or periodic BC on the sides.
  int i, j, k, im, ip, jm, jp, km, kp, nxy, newnxy;
  int newindex;
  int indpmm, indpmp, indppm, indppp;
  int indmmm, indmmp, indmpm, indmpp;
  nxy = phi->nx * phi->ny;
  newnxy = newphi->nx * newphi->ny;
  for (i=0; i<newphi->nx; i++)
    {
      im = max(0,rint((float)i / 2.0 - 0.1));
      ip = min(newphi->nx-1,rint((float)i / 2.0 + 0.1));
      for (j=0; j<newphi->ny; j++)
	{
	  jm = max(0,rint((float)j / 2.0 - 0.1));
	  jp = min(newphi->ny-1,rint((float)j / 2.0 + 0.1));
	  if (newBCType->data[i + j * newphi->nx] == 1) // Free BC at z = 0
	    {
	      newindex = i + j * newphi->nx;
	      indmmm = im + jm * phi->nx;
	      indmmp = im + jm * phi->nx;
	      indmpm = im + jp * phi->nx;
	      indmpp = im + jp * phi->nx;
	      indpmm = ip + jm * phi->nx;
	      indpmp = ip + jm * phi->nx;
	      indppm = ip + jp * phi->nx;
	      indppp = ip + jp * phi->nx;
	      newphi->data[newindex] = (phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp]) / 8.0;
	    }
	  for (k=newphi->Vkmin[i]; k<newphi->nz-1; k++)
	    {
	      km = rint((float)k / 2.0 - 0.1);
	      kp = rint((float)k / 2.0 + 0.1);
	      newindex = i + j * newphi->nx + k * newnxy;
	      indmmm = im + jm * phi->nx + km * nxy;
	      indmmp = im + jm * phi->nx + kp * nxy;
	      indmpm = im + jp * phi->nx + km * nxy;
	      indmpp = im + jp * phi->nx + kp * nxy;
	      indpmm = ip + jm * phi->nx + km * nxy;
	      indpmp = ip + jm * phi->nx + kp * nxy;
	      indppm = ip + jp * phi->nx + km * nxy;
	      indppp = ip + jp * phi->nx + kp * nxy;
	      newphi->data[newindex] = (phi->data[indmmm] + phi->data[indmmp] + phi->data[indmpm] + phi->data[indmpp] + phi->data[indpmm] + phi->data[indpmp] + phi->data[indppm] + phi->data[indppp]) / 8.0;
	    }
	}
    }
  return;
}



void MultiGrid::VCycle_Inner(Array3D** phi, Array3D** rho, Array3D** elec, Array3D** hole, Array3D** eps, Array2D** BCType,Array2D** QFe, Array2D** QFh, double w, int nsteps, int ncycle)
{
  // iterates the coarsest grid down to machine precision,
  // then does a given number of steps (ncycle) at each finer scale on the way up.
  int i, j, niter, NumSOR;
  double error, AveIterations;
  for (i=0; i<nsteps; i++)
    {
      Restrict(phi[i], phi[i+1], rho[i], rho[i+1], BCType[i], BCType[i+1]);
    }
  //CountCharges(rho, elec, hole);
  error = 100.0; niter = 0;
  NumSOR = 0;
  AveIterations = 0;
  while (error > 1.0E-12)
    {
      niter += 1;
      AveIterations = SOR_Inner(phi[nsteps], rho[nsteps], elec[nsteps], hole[nsteps], eps[nsteps], BCType[nsteps], QFe[nsteps], QFh[nsteps], w);
      NumSOR += 1;
      error = Error_Inner(phi[nsteps], rho[nsteps], elec[nsteps], hole[nsteps], eps[nsteps]);      
    }
  if(VerboseLevel > 1)
    {
      printf("Completed iterations at resolution %dx%dx%d. Number of steps = %d. Error = %.3g\n",phi[nsteps]->nx-1,phi[nsteps]->ny-1,phi[nsteps]->nz-1,niter,error);
      fflush(stdout);
    }
  for (i=nsteps; i>0; i--)
    {
      Prolongate(phi[i], phi[i-1], BCType[i-1]);
      if (i == 1) niter = ncycle;
      else niter = ncycle * 8 * (int)pow(2,i-1);
      NumSOR = 0;
      AveIterations = 0;
      for (j=0; j<niter; j++)
	{
	  AveIterations = SOR_Inner(phi[i-1], rho[i-1], elec[i-1], hole[i-1], eps[i-1], BCType[i-1], QFe[i-1], QFh[i-1], w);
	  NumSOR += 1;
	}
      error = Error_Inner(phi[i-1], rho[i-1], elec[i-1], hole[i-1], eps[i-1]);      
      if(VerboseLevel > 1)
	{
	  printf("Completed iterations at resolution %dx%dx%d. Number of steps = %d. Error = %.3g Ave inner iterations = %.2f\n",phi[i-1]->nx-1,phi[i-1]->ny-1,phi[i-1]->nz-1,j,error,AveIterations/(double)NumSOR);
	  fflush(stdout);
	}
    }
  return;
}

void MultiGrid::VCycle_Zero(Array3D** phi, Array3D** rho, Array3D** elec, Array3D** hole, Array3D** eps, Array2D** BCType,Array2D** QFe, Array2D** QFh, double w, int ncycle)
{
  // Only does a given number of steps (ncycle) at the finest scale
  int j, NumSOR;
  double error, AveIterations;
  NumSOR = 0;
  AveIterations = 0;
      for (j=0; j<ncycle; j++)
    {
      AveIterations = SOR_Inner(phi[0], rho[0], elec[0], hole[0], eps[0], BCType[0], QFe[0], QFh[0], w);
      NumSOR += 1;
    }
  error = Error_Inner(phi[0], rho[0], elec[0], hole[0], eps[0]);  
  if(VerboseLevel > 1)
    {
      printf("Completed iterations at resolution %dx%dx%d. Number of steps = %d. Error = %.3g Ave inner iterations = %.2f\n",phi[0]->nx-1,phi[0]->ny-1,phi[0]->nz-1,j,error,AveIterations/(double)NumSOR);
      fflush(stdout);
    }
  return;
}

void MultiGrid::WriteOutputFile(string outputfiledir, string filenamebase, string name, Array3D* array)
{
  // This writes the data to the HDF files
  string underscore = "_", slash = "/", hdfname, filename;
  int* int_attr_data  = new int[3];
  double* double_attr_data  = new double[3];
  double* flipped_data  = new double[array->nx*array->ny*array->nz];
  int i, j, k, index, flipped_index;

  // There must be a better way to change from x fast to z fast, but this works.
  for (i=0; i<array->nx; i++)
    {
      for (j=0; j<array->ny; j++)
	{
	  for (k=0; k<array->nz; k++)
	    {
	      index = i + j * array->nx + k * array->nx * array->ny;
	      flipped_index = k + j * array->nz + i * array->nz * array->ny;
	      flipped_data[flipped_index] = array->data[index];
	    }
	}
    }
  hdfname = filenamebase+underscore+name;
  filename = outputfiledir+slash+hdfname+".hdf5";
  WriteHDF5File3(filename, hdfname, array->nx, array->ny, array->nz, flipped_data);
  // Now we write the attributes
  int_attr_data[0] = array->nx; int_attr_data[1] = array->ny; int_attr_data[2] = array->nz;
  WriteHDF5IntAttribute(filename, hdfname, "Dimension", 3, int_attr_data);
  double_attr_data[0] = array->xmin; double_attr_data[1] = array->ymin; double_attr_data[2] = array->zmin;
  WriteHDF5DoubleAttribute(filename, hdfname, "Lower_Left", 3, double_attr_data);
  double_attr_data[0] = array->xmax; double_attr_data[1] = array->ymax; double_attr_data[2] = array->zmax;
  WriteHDF5DoubleAttribute(filename, hdfname, "Upper_Right", 3, double_attr_data);

  delete[] flipped_data;
  return;
}

void MultiGrid::ReadOutputFile(string outputfiledir, string filenamebase, string name, Array3D* array)
{
  // This reads the data from the HDF files
  string underscore = "_", slash = "/", hdfname, filename;
  double* flipped_data  = new double[array->nx*array->ny*array->nz];
  int i, j, k, index, flipped_index;
  hdfname = filenamebase+underscore+name;
  filename = outputfiledir+slash+hdfname+".hdf5";
  ReadHDF5File3(filename, hdfname, array->nz, array->ny, array->nx, flipped_data);
  // There must be a better way to change from x fast to z fast, but this works.
  for (i=0; i<array->nx; i++)
    {
      for (j=0; j<array->ny; j++)
	{
	  for (k=0; k<array->nx; k++)
	    {
	      index = i + j * array->nx + k * array->nx * array->ny;
	      flipped_index = k + j * array->nz + i * array->nz * array->ny;
	      array->data[index] = flipped_data[flipped_index];
	    }
	}
    }
  delete[] flipped_data;
  return;
}

void MultiGrid::Gradient(Array3D* phi, Array3D** E)
{
  int i, j, k, nxy, ind;
  nxy = phi->nx * phi->ny;

  // Ex
  for (j=0; j<phi->ny; j++)
    {
      for (k=0; k<phi->nz; k++)
	{
	  ind = j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind]) / phi->dx;
	  ind = 1 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind-1]) / (2.0 * phi->dx);
	  ind = phi->nx-2 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind+1] - phi->data[ind-1]) / (2.0 * phi->dx);
	  ind = phi->nx-1 + j * phi->nx + k * nxy;
	  E[0]->data[ind] = (phi->data[ind] - phi->data[ind-1]) / phi->dx;
	  for (i=2; i<phi->nx-2; i++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[0]->data[ind] = (-phi->data[ind+2] + 8.0 * phi->data[ind+1] - 8.0 * phi->data[ind-1] + phi->data[ind-2]) / (12.0 * phi->dx);
	    }
	}
    }

    // Ey
  for (i=0; i<phi->nx; i++)
    {
      for (k=0; k<phi->nz; k++)
	{
	  ind = i + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind]) / phi->dy;
	  ind = i + phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind-phi->nx]) / (2.0 * phi->dy);
	  ind = i + (phi->ny-2) * phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind+phi->nx] - phi->data[ind-phi->nx]) / (2.0 * phi->dy);
	  ind = i + (phi->ny-1) * phi->nx + k * nxy;
	  E[1]->data[ind] = (phi->data[ind] - phi->data[ind-phi->nx]) / phi->dy;
	  for (j=2; j<phi->ny-2; j++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[1]->data[ind] = (-phi->data[ind+2*phi->nx] + 8.0 * phi->data[ind+phi->nx] - 8.0 * phi->data[ind-phi->nx] + phi->data[ind-2*phi->nx]) / (12.0 * phi->dy);
	    }
	}
    }

  // Ez
  for (i=0; i<phi->nx; i++)
    {
      for (j=0; j<phi->ny; j++)
	{
	  ind = i + j * phi->nx;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind]) / phi->dzp * phi->dzpdz[0];
	  ind = i + j * phi->nx + nxy;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind-nxy]) / (2.0 * phi->dzp)  * phi->dzpdz[1];
	  ind = i + j * phi->nx + (phi->nz-2) * nxy;
	  E[2]->data[ind] = (phi->data[ind+nxy] - phi->data[ind-nxy]) / (2.0 * phi->dzp) * phi->dzpdz[phi->nz-2];
	  ind = i + j * phi->nx + (phi->nz-1) * nxy;
	  E[2]->data[ind] = (phi->data[ind] - phi->data[ind-nxy]) / phi->dzp * phi->dzpdz[phi->nz-1];
	  for (k=2; k<phi->nz-2; k++)
	    {
	      ind = i + j * phi->nx + k * nxy;
	      E[2]->data[ind] = (-phi->data[ind+2*nxy] + 8.0 * phi->data[ind+nxy] - 8.0 * phi->data[ind-nxy] + phi->data[ind-2*nxy]) / (12.0 * phi->dzp) * phi->dzpdz[k];
	    }
	}
    }
  return;
}

void MultiGrid::Trace(double* point, int bottomsteps, bool savecharge, double bottomcharge, ofstream& file)
{
  // This traces an electron down to the bottom, saving path info if requested
  // Diffusion has now been added. This version recalculates mu at each point.
  // And iterates bottomsteps steps after reaching the bottom.
  // If savecharge is true, it finds and stores the self-consistent charge locations
  // If ElectronAccumulation=0, it will save bottomcharge in each location for bottomsteps steps.
  // If ElectronAccumlate=0, it will only determine the pixel and relay on QFe to determine the charge location.
  // id is a unique integer identifier for each electron that we track.
  // phase encodes the tracking phase and is recorded to the Pts file when LogPixelPaths != 0.
  // 0 - initial position.
  // 1 - endpoint of a step through the bulk.
  // 2 - just reached bottom and settling to equilibrium.
  // 3 - equilibrium motion at the bottom, logging charge.
  // 4 - final position.

  int i, j, k, nsteps = 0, nstepsmax = 10000;
  bool ReachedBottom = false;
  double mu, E2, Emag, ve, vth, tau, Tscatt;
  double theta, phiangle, zmin, zbottom;
  double x, y, z;

  zmin = E[0]->z[Channelkmax];
  zbottom = E[0]->zmz[Channelkmin] + 0.01;
  x = point[0]; y = point[1]; z = point[2];
  double*  E_interp = new double[3];
  E2 = 0.0;
  for (i=0; i<3; i++)
    {
      E_interp[i] = E[i]->DataInterpolate3D(point[0],point[1],point[2]);
      E2 += E_interp[i] * E_interp[i];
    }
  Emag = max(0.1, sqrt(E2));
  mu = mu_Si(Emag * MICRON_PER_CM, CCDTemperature); // Mobility
  vth = sqrt(3.0 * KBOLTZMANN * CCDTemperature / ME)  * MICRON_PER_M * DiffMultiplier; // Thermal Velocity
  vth = vth / sqrt((double)NumDiffSteps);
  tau  = ME / QE * mu * METER_PER_CM * METER_PER_CM; // scattering time

  static int id = 0;
  int phase = 0;
  // Log initial position.
  file << setw(8) << id << setw(8) << nsteps << setw(3) << phase
       << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;

  phase = 1;
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
      phiangle = 2.0 * pi * drand48();
      theta = acos(-1.0 + 2.0 * drand48());
      Tscatt = -tau * log(1.0 - drand48()) * (double)NumDiffSteps;
      point[0] += (vth * sin(theta) * cos(phiangle) + E_interp[0] * ve) * Tscatt;
      point[1] += (vth * sin(theta) * sin(phiangle) + E_interp[1] * ve) * Tscatt;
      point[2] += (vth * cos(theta) + E_interp[2] * ve) * Tscatt;

      if (point[2] < zmin && !ReachedBottom)
	{
	  ReachedBottom = true;
	  nstepsmax = nsteps + bottomsteps + EquilibrateSteps;
	  phase = 2;
	  // After reaching bottom, iterate bottomsteps (*1.1) more steps.
	  // The first bottomsteps/10 steps are to let it settle to an
	  // equilibrium location, then we start logging the charge location
	}
      if (ReachedBottom && nsteps > nstepsmax - bottomsteps)
	{
	  //  Start logging location after EquilibrateSteps.
	  phase = (nsteps < nstepsmax) ? 3 : 4;
	  if (point[2] <= zbottom && SaturationModel == 1)
	    {
	      break; // Electron recombines and is lost if it reaches the gate oxide
	    }
	  point[2] = max(zbottom, point[2]);
	  i = E[0]->XIndex(point[0]);
	  j = E[0]->YIndex(point[1]);
	  k = E[0]->ZIndex(point[2]);
	  if ((hole[0]->data[i + j * hole[0]->nx + k * hole[0]->nx * hole[0]->ny] > 0.1) && SaturationModel == 1)
	    {
	      break; // Electron recombines if it encounters free holes.
	    }
	  if (savecharge && ElectronAccumulation == 1)
	    {
	      // Find the pixel the charge is in and add 1 electron to it.
	      int PixX = (int)floor((point[0] - PixelBoundaryLowerLeft[0]) / PixelSize);
	      int PixY = (int)floor((point[1] - PixelBoundaryLowerLeft[1]) / PixelSize);
	      j = PixX + PixelBoundaryNx * PixY;
	      if (j > 0 && j < PixelBoundaryNx * PixelBoundaryNy)   CollectedCharge[0][j] += 1;
	      break;
	    }
	  if (i > 0 && i < elec[0]->nx-1 && j > 0 && j < elec[0]->ny-1 && k < elec[0]->nz-1 && savecharge && ElectronAccumulation == 0)
	    {
	      elec[0]->data[i + j * elec[0]->nx + k * elec[0]->nx * elec[0]->ny] += bottomcharge;// Add bottomcharge to this grid cell
	    }
	}
      if(LogPixelPaths == 1) 
      {
        // Log latest position update.
        file << setw(8) << id << setw(8) << nsteps << setw(3) << phase
	     << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
      }
    point[2] = max(zbottom, point[2]);
    }
  delete[] E_interp;
  
  if(LogPixelPaths == 0) 
    {
      // Log final position
      file << setw(8) << id << setw(8) << nsteps << setw(3) << phase
	   << setw(15) << point[0] << setw(15) << point[1] << setw(15) << point[2] << endl;
    }
  id += 1;
  return;
}


double MultiGrid::GetElectronInitialZ() {
    if(FilterIndex >= 0 && FilterIndex < n_band) {
        int cdf_index = (int)floor(n_filter_cdf * drand48());
        return SensorThickness - filter_cdf[FilterIndex * n_filter_cdf + cdf_index];
    }
    else {
        return ElectronZ0Fill;
    }
}

void MultiGrid::TraceSpot(int m)
{
  // This builds up a Gaussian spot with given center (Xoffset, Yoffset) and SigmaX and SigmaY
  double x, y, z, rsq, v1, v2, fac, xcenter, ycenter;
  int n;
  double bottomcharge = 1.0 / (double)BottomSteps;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = boost::lexical_cast<std::string>(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;

  xcenter = (PixelBoundaryUpperRight[0] + PixelBoundaryLowerLeft[0]) / 2.0 + Xoffset;
  ycenter = (PixelBoundaryUpperRight[1] + PixelBoundaryLowerLeft[1]) / 2.0 + Yoffset;
  for (n=0; n<NumElec; n++)
    {
      //  Use Box-Muller algorithm to generate two Gaussian random numbers
      rsq = 1000.0;
      while (rsq >= 1.0 || rsq == 0.0)
	{
	  v1 = 2.0 * drand48() - 1.0;
	  v2 = 2.0 * drand48() - 1.0;
	  rsq = v1*v1 + v2 *v2;
	}
      fac = sqrt(-2.0 * log(rsq) / rsq);
      x = xcenter + Sigmax * v1 * fac;
      y = ycenter + Sigmay * v2 * fac;
      point[0] = x;
      point[1] = y;
      z = GetElectronInitialZ();
      point[2] = z;
      Trace(point, BottomSteps, true, bottomcharge, file);
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;
  return;
}

void MultiGrid::TraceList(int m)
{
  // This builds up a Gaussian spot with given center (Xoffset, Yoffset) and SigmaX and SigmaY
  double x, y, z, zbottom, xcenter, ycenter, abs_length, path_length;
  zbottom = E[0]->zmz[Channelkmin] + 0.01;
  int n, nlist;
  double bottomcharge = 1.0 / (double)BottomSteps;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = boost::lexical_cast<std::string>(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;

  xcenter = (PixelBoundaryUpperRight[0] + PixelBoundaryLowerLeft[0]) / 2.0 + Xoffset;
  ycenter = (PixelBoundaryUpperRight[1] + PixelBoundaryLowerLeft[1]) / 2.0 + Yoffset;
  for (n=0; n<NumElec; n++)
    {
      nlist = m * NumElec + n;
      if (nlist > NumPhotons)
	{
	  printf("Reached end of photon list in step %d.\n",m);
	  return;
	}

      abs_length = pow(10.0,(-4.0 + (PhotonListlambda[nlist] - 500.0) / 250.0)) * 1.0E4; //Approximate formula in micron^-1
      path_length = -abs_length * log(1.0 - drand48());
      x = xcenter + PixelSize * PhotonListx[nlist]; // in microns
      y = ycenter + PixelSize * PhotonListy[nlist]; // in microns      
      x += PhotonListdxdz[nlist] * path_length;
      y += PhotonListdydz[nlist] * path_length;
      z = SensorThickness - path_length / sqrt(1.0 + PhotonListdxdz[nlist]*PhotonListdxdz[nlist] + PhotonListdydz[nlist]*PhotonListdydz[nlist]); // in microns
      if (x < PixelBoundaryLowerLeft[0] || x > PixelBoundaryUpperRight[0] || y < PixelBoundaryLowerLeft[1] || y > PixelBoundaryUpperRight[1] || z > SensorThickness || z < zbottom) continue;
      if (nlist % 1000 == 0)
	{
	  printf("Nlist = %d, (x,y,z) = (%f,%f,%f)\n",nlist,x,y,z);
	}
      point[0] = x;
      point[1] = y;
      point[2] = z;
      Trace(point, BottomSteps, true, bottomcharge, file);
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;
  return;
}

void MultiGrid::TraceGrid(int m)
{
  // This traces a grid of starting electron locations.
  double x, y, z;
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = boost::lexical_cast<std::string>(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;

  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  x = PixelBoundaryLowerLeft[0] + PixelBoundaryStepSize[0] / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      y = PixelBoundaryLowerLeft[1] + PixelBoundaryStepSize[1] / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  point[0] = x;
	  point[1] = y;
	  z = GetElectronInitialZ();
	  point[2] = z;
	  Trace(point, 100, false, 0.0, file);
	  y += PixelBoundaryStepSize[1];
	}
      x += PixelBoundaryStepSize[0];
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;
  return;
}

void MultiGrid::TraceRegion(int m)
{
  // This traces a random set of starting electron locations within the PixelBoundary.
  double x, y, z, boxx, boxy;
  int n;
  boxx = PixelBoundaryUpperRight[0] - PixelBoundaryLowerLeft[0];
  boxy = PixelBoundaryUpperRight[1] - PixelBoundaryLowerLeft[1];
  double* point = new double[3];
  string underscore = "_", slash = "/", name = "Pts";
  string StepNum = boost::lexical_cast<std::string>(m);
  string filename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+name+".dat");
  ofstream file;
  file.open(filename.c_str());
  file.setf(ios::fixed);
  file.setf(ios::showpoint);
  file.setf(ios::left);
  file.precision(4);
  // Write header line.
  file << setw(8) << "id" << setw(8) << "step" << setw(3) << "ph"
       << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
  // If LogPixelPaths > 1, then we only log some pixel paths
  int OldLogPixelPaths = LogPixelPaths;
  if (LogPixelPaths != 0)
    {
      if (m % LogPixelPaths == 0)
	{
	  LogPixelPaths = 1;
	}
      else
	{
	  LogPixelPaths = 0;
	}
    }
  // Initialize for PixelBoundaryTestType == 4 mode.
  double bottomcharge = 1.0 / (double)BottomSteps;
  double x_center = 0.5 * (PixelBoundaryLowerLeft[0] + PixelBoundaryUpperRight[0]);
  double y_center = 0.5 * (PixelBoundaryLowerLeft[1] + PixelBoundaryUpperRight[1]);

  for (n=0; n<NumElec; n++)
    {
      if (n%1000==0)
	{
	  printf("Finished %d electrons\n",n);
	  fflush(stdout);
	}
      if(PixelBoundaryTestType == 4) {
          x = x_center + (drand48() - 0.5) * PixelSize;
          y = y_center + (drand48() - 0.5) * PixelSize;
          z = GetElectronInitialZ();
      }
      else {
          x = PixelBoundaryLowerLeft[0] + drand48() * boxx;
          y = PixelBoundaryLowerLeft[1] + drand48() * boxy;
          z = ElectronZ0Fill;
      }
      point[0] = x;
      point[1] = y;
      point[2] = z;
      if(PixelBoundaryTestType == 4) {
          // Accumulate charge, the same as TraceSpot().
          Trace(point, BottomSteps, true, bottomcharge, file);
      }
      else {
          // Do not accumulate charge, for backwards compatibility.
          Trace(point, 100, false, 0.0, file);
      }
    }
  file.close();
  printf("Finished writing grid file - %s\n",filename.c_str());
  fflush(stdout);
  delete[] point;
  LogPixelPaths = OldLogPixelPaths;  
  return;
}

void MultiGrid::FindEdge(double* point, double theta, ofstream& file)
{
  // This finds the edge of the pixel through binary search given a starting point and a line angle
  int nsteps, pixx, pixy, newpixx, newpixy;
  double sinth, costh, x, y, z0, deltar, tolerance;
  sinth = sin(theta);
  costh = cos(theta);
  z0 = point[2];
  x = point[0]; y = point[1];
  pixx = (int)floor((point[0] - PixelBoundaryLowerLeft[0]) / PixelSize);
  pixy = (int)floor((point[1] - PixelBoundaryLowerLeft[1]) / PixelSize);
  deltar = 0.5 * PixelSize;
  tolerance = 0.0001;
  nsteps = 0;
  while (fabs(deltar) > tolerance)
    {
      nsteps += 1;
      if (nsteps >200)
	{
	  printf("Too many steps in edge finding, pixx = %d, pixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",pixx, pixy, theta, nsteps,x,y);
	  fflush(stdout);
	  break;
	}
      x += deltar * costh;
      y += deltar * sinth;
      point[0] = x;
      point[1] = y;
      point[2] = z0;
      Trace(point, 100, false, 0.0, file);
      newpixx = (int)floor((point[0] - PixelBoundaryLowerLeft[0]) / PixelSize);
      newpixy = (int)floor((point[1] - PixelBoundaryLowerLeft[1]) / PixelSize);
      //printf("Finding edge, newpixx = %d, newpixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",newpixx, newpixy, theta, nsteps,point[0],point[1]);
      if (newpixx != pixx || newpixy != pixy)
	{
	  x -= deltar * costh;
	  y -= deltar * sinth;
	  deltar /= 2.0;
	}
    }
  //printf("Found edge, pixx = %d, pixy = %d, theta = %.3f, %d steps, x = %.3f, y = %.3f\n",pixx, pixy, theta, nsteps,x,y);
  fflush(stdout);
  point[0] = x;
  point[1] = y;
  point[2] = z0;
  return;
}

void MultiGrid::FindCorner(double* point, double* theta, ofstream& file)
{
  // This finds the corner of the pixel through binary search given a starting point and a line angle
  int nsteps;
  double theta0, delta_theta, r, rm, r0, rp, deltar, tolerance, x0, y0, z0;
  theta0 = *theta;
  delta_theta = 0.10;
  z0 = point[2];
  x0 = point[0]; y0 = point[1];
  deltar = 1.0;
  tolerance = 0.0001;
  nsteps = 0;
  while (fabs(deltar) > tolerance)
    {
      nsteps += 1;
      if (nsteps >200)
	{
	  printf("Too many steps in corner finding, theta = %.3f, %d steps, x0 = %.3f, y0 = %.3f\n",theta0, nsteps,x0,y0);
	  fflush(stdout);
	  break;
	}
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      r0 = r;
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0 + delta_theta, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      rp = r;
      point[0] = x0; point[1] = y0; point[2] = z0;
      FindEdge(point, theta0 - delta_theta, file);
      r = sqrt((point[0] - x0) * (point[0] - x0) + (point[1] - y0) * (point[1] - y0));
      rm = r;
      //printf("Step = %d, theta = %.5f, r0 = %.4f, rm = %.4f, rp = %.4f\n",nsteps,theta,r0,rm,rp);
      if (rp > r0)
	{
	  theta0 = theta0 + delta_theta;
	  deltar = fabs(r0 - rp);
	}
      else if (rm > r0)
	{
	  theta0 = theta0 - delta_theta;
	  deltar = fabs(r0 - rm);
	}
      else if (rp > rm)
	{
	  delta_theta = delta_theta / 2.0;
	  deltar = fabs(r0 - rp);
	}
      else
	{
	  delta_theta = delta_theta / 2.0;
	  deltar = fabs(r0 - rm);
	}
    }
  //printf("Found corner, theta = %.4f, x = %.4f, y = %.4f, after %d steps\n",theta0, point[0], point[1], nsteps);
  //fflush(stdout);
  *theta = theta0;
  return ;
}

void MultiGrid::CalculatePixelAreas(int m)
{
  // This finds the pixel vertices and the areas of the pixel grid.
  int OldLogPixelPaths, k, n, pixx, pixy;
  //Turn off Pixel Path logging
  OldLogPixelPaths = LogPixelPaths;
  LogPixelPaths = 0;
  //Turn off Diffusion
  double OldDiffMultiplier;
  OldDiffMultiplier = DiffMultiplier;
  DiffMultiplier = 0.0;
  double x, y, xb, yb, theta, theta0, theta1, dtheta, area;
  double* point = new double[3];
  string dummy = "dummy";
  string ptsfilename = (dummy);
  ofstream ptsfile; //Not needed, but we need to pass something to the Trace subroutine


  // Now calculate the pixel vertices
  Polygon** polyarray = new Polygon*[PixelBoundaryNx * PixelBoundaryNy];
  x = PixelBoundaryLowerLeft[0] + PixelSize / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      pixx = (int)floor((x - PixelBoundaryLowerLeft[0]) / PixelSize);
      y = PixelBoundaryLowerLeft[1] + PixelSize / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  pixy = (int)floor((y - PixelBoundaryLowerLeft[1]) / PixelSize);
	  polyarray[pixx + PixelBoundaryNx * pixy] = new Polygon(4 * NumVertices + 4);
	  //First, find the four corners
	  for (n=1; n<8; n+=2)
	    {
	      theta = (double)n * pi / 4.0;
	      point[0] = x;
	      point[1] = y;
	      point[2] = ElectronZ0Area;
	      FindCorner(point, &theta, ptsfile);
	      Point* two_d_point = new Point(point[0], point[1], theta);
	      polyarray[pixx + PixelBoundaryNx * pixy]->AddPoint(two_d_point);
	    }
	  // Now find NumVertices points along each edge
	  for (n=0; n<4; n++) // Four corners
	    {
	      for (k=0; k<NumVertices; k++)
		{
		  if (n == 3) // Need this to make angles continuous through zero
		    {
		      theta0 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta - 2.0 * pi;
		      theta1 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->theta;
		    }
		  else
		    {
		      theta0 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta;
		      theta1 = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[(n+1)]->theta;
		    }
		  dtheta = (theta1 - theta0) / ((double)NumVertices + 1.0);
		  theta  = theta0  + ((double)k + 1.0) * dtheta;
		  point[0] = x;
		  point[1] = y;
		  point[2] = ElectronZ0Area;
		  FindEdge(point, theta, ptsfile);
		  //printf("Found edge, pixx = %d, pixy = %d, theta = %.3f, x = %.3f, y = %.3f\n",pixx, pixy, theta, point[0],point[1]);
		  Point* two_d_point = new Point(point[0], point[1], theta);
		  polyarray[pixx + PixelBoundaryNx * pixy]->AddPoint(two_d_point);
		}
	    }
	  printf("Finished vertex finding for pixel, pixx = %d, pixy = %d\n",pixx, pixy);
	  printf("Corners at (%.3f, %.3f),(%.3f, %.3f),(%.3f, %.3f),(%.3f, %.3f)\n",
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[0]->y,
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[1]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[1]->y,
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[2]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[2]->y,
		 polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[3]->x, polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[3]->y);
	  fflush(stdout);
	  y += PixelSize;
	}
      x += PixelSize;
    }

  // Now calculate and print out the pixel areas
  string underscore = "_", slash = "/", vertexname = "Vertices", areaname = "Area";
  string StepNum = boost::lexical_cast<std::string>(m);
  string areafilename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+areaname+".dat");
  ofstream areafile;
  areafile.open(areafilename.c_str());
  areafile.setf(ios::fixed);
  areafile.setf(ios::showpoint);
  areafile.setf(ios::left);
  areafile.precision(4);
  areafile  << setw(15) << "Nx" << setw(15) << "Ny" << setw(15) << "Area" << endl;
  for (pixx=0; pixx<PixelBoundaryNx; pixx++)
    {
      for (pixy=0; pixy<PixelBoundaryNy; pixy++)
	{
	  area = polyarray[pixx + PixelBoundaryNx * pixy]->Area();
	  printf("Found area, pixx = %d, pixy = %d, area = %.3f\n",pixx, pixy, area);
	  fflush(stdout);
	  areafile  << setw(15) << pixx << setw(15) << pixy << setw(15) << area << endl;
	}
    }
  areafile.close();
  printf("Finished writing grid file - %s\n",areafilename.c_str());

  // Now print out the pixel vertices
  string vertexfilename = (outputfiledir+slash+outputfilebase+underscore+StepNum+underscore+vertexname+".dat");
  ofstream vertexfile;
  vertexfile.open(vertexfilename.c_str());
  vertexfile.setf(ios::fixed);
  vertexfile.setf(ios::showpoint);
  vertexfile.setf(ios::left);
  vertexfile.precision(4);
  vertexfile  << setw(15) << "X0" << setw(15) << "Y0" << setw(15)<< "Theta" << setw(15) << "X" << setw(15) << "Y" << endl;
  x = PixelBoundaryLowerLeft[0] + PixelSize / 2.0;
  while (x < PixelBoundaryUpperRight[0])
    {
      pixx = (int)floor((x - PixelBoundaryLowerLeft[0]) / PixelSize);
      y = PixelBoundaryLowerLeft[1] + PixelSize / 2.0;
      while (y < PixelBoundaryUpperRight[1])
	{
	  pixy = (int)floor((y - PixelBoundaryLowerLeft[1]) / PixelSize);
	  for (n=0; n<polyarray[pixx + PixelBoundaryNx * pixy]->npoints; n++)
	    {
	      xb = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->x;
	      yb = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->y;
	      theta = polyarray[pixx + PixelBoundaryNx * pixy]->pointlist[n]->theta;
	      vertexfile  << setw(15) << x << setw(15) << y << setw(15)<< theta << setw(15)<< xb << setw(15) << yb << endl;
	    }
	  y += PixelSize;
	}
      x += PixelSize;
    }
  vertexfile.close();
  ptsfile.close();
  printf("Finished writing grid file - %s\n",vertexfilename.c_str());
  fflush(stdout);
  // Clean up
  delete[] point;
  for (pixx=0; pixx<PixelBoundaryNx; pixx++)
    {
      for (pixy=0; pixy<PixelBoundaryNy; pixy++)
	{
	  delete polyarray[pixx + PixelBoundaryNx * pixy];
	}
    }
  delete[] polyarray;
  //Put Pixel path logging back where it was
  LogPixelPaths = OldLogPixelPaths;
  //Turn Diffusion back on
  DiffMultiplier = OldDiffMultiplier;
  return;
}

double MultiGrid::mu_Si (double E,double T)
{
  // Shamelssly copied from phosim
  // Jacobini et al. (1977) equation (9)
  double vm=1.53e9 * pow(T,-0.87); // cm/s
  double Ec = 1.01 * pow(T,1.55); // V/cm
  double beta = 2.57e-2 * pow(T,0.66); // index
  return((vm/Ec)/pow(1 + pow(fabs(E)/Ec,beta),1/beta));
}

void MultiGrid::Set_QFh(Array2D** QFh, double value)
{
  //Set hole quasi-Fermi level in the region to give the desired number of electrons
  int i, j, n, index;
  for (n=0; n<nsteps+1; n++)
    {
      for (i=0; i<QFh[n]->nx; i++)
	{
	  for (j=0; j<QFh[n]->ny; j++)
	    {
	      index = i + j * QFh[n]->nx;
	      QFh[n]->data[index] = value;
	    }
	}
    }
  printf("QFh = %.2f\n",value);
  fflush(stdout);
  return;
  }

void MultiGrid::Adjust_QFe(Array2D** QFe, Array3D* rho, Array3D* elec)
{
  // Set electron quasi-Fermi level in the region to give the desired number of electrons
  int i, j, k, n, m, q, index, index2, imin, imax, jmin, jmax, PixX, PixY;
  double PixXmin, PixXmax, PixYmin, PixYmax, TotalElectrons = 0.0;
  // First clear the array
  for (n=0; n<nsteps+1; n++)
    {
      for (i=0; i<QFe[n]->nx; i++)
	{
	  for (j=0; j<QFe[n]->ny; j++)
	    {
	      index = i + j * QFe[n]->nx;
	      QFe[n]->data[index] = 40.0;
	    }
	}
    }

  // Now count the electrons in the pixel regions.
  for (m=0; m<NumberofPixelRegions; m++)
    {
      for (q=0; q<NumberofFilledWells[m]; q++)
	{
	  // First figure out the pixel coordinates
	  PixX = (int)floor((FilledPixelCoords[m][q][0] - PixelRegionLowerLeft[m][0]) / PixelSize);
	  PixY = (int)floor((FilledPixelCoords[m][q][1] - PixelRegionLowerLeft[m][1]) / PixelSize);
	  PixXmin = PixelRegionLowerLeft[m][0] + (double)PixX * PixelSize;
	  PixYmin = PixelRegionLowerLeft[m][1] + (double)PixY * PixelSize;
	  PixXmax = PixXmin + PixelSize;
	  PixYmax = PixYmin + PixelSize;	      
	  imin = elec->XIndex(PixXmin);
	  imax = elec->XIndex(PixXmax);
	  jmin = elec->YIndex(PixYmin);
	  jmax = elec->YIndex(PixYmax);
	  // Then count how many electrons currently in the pixel
	  ElectronCount[m][q] = 0.0;
	  for (i=imin; i<imax; i++)
	    {
	      for (j=jmin; j<jmax; j++)
		{
		  for (k=rho->Ckmin[i]; k<elec->nz; k++)
		    {
		      index = i + j * elec->nx + k * elec->nx * elec->ny;
		      ElectronCount[m][q] += elec->data[index];
		      TotalElectrons +=  elec->data[index];
		      if (BuildQFeLookup == 1)
			{
			  elec->data[index] = 0.0; // Clear the electrons if building the lookup table
			}
		    }
		}
	    }
	  //printf("m=%d,q=%d,PixX = %.1f, Pixy = %.1f, Desired Electrons = %d, ElecCount = %.4f, QFe = %.4f\n",m,q,FilledPixelCoords[m][q][0], FilledPixelCoords[m][q][1], CollectedCharge[m][q], ElectronCount[m][q], PixelQFe[m][q]);
	  if (BuildQFeLookup == 1)
	    {
	      PixelQFe[m][q] = QFemin + (QFemax - QFemin) / (double)NQFe * (double)q;
	      QFeLookup[q] = ElectronCount[m][q];
	      //printf("PixX = %.1f, Pixy = %.1f, ElecCount = %.4f, QFe = %.4f\n",FilledPixelCoords[m][q][0], FilledPixelCoords[m][q][1], ElectronCount[m][q], PixelQFe[m][q]);
	    }
	  else
	    {
	      PixelQFe[m][q] = ElectronQF(CollectedCharge[m][q]);
	    }
	  // Now set pixel QFe to calculated value for all steps
	  for (n=0; n<nsteps+1; n++)
	    {
	      // First figure out the pixel coordinates
	      PixX = (int)floor((FilledPixelCoords[m][q][0] - PixelRegionLowerLeft[m][0]) / PixelSize);
	      PixY = (int)floor((FilledPixelCoords[m][q][1] - PixelRegionLowerLeft[m][1]) / PixelSize);
	      PixXmin = PixelRegionLowerLeft[m][0] + (double)PixX * PixelSize;
	      PixYmin = PixelRegionLowerLeft[m][1] + (double)PixY * PixelSize;
	      PixXmax = PixXmin + PixelSize;
	      PixYmax = PixYmin + PixelSize;	      
	      imin = QFe[n]->XIndex(PixXmin);
	      imax = QFe[n]->XIndex(PixXmax);
	      jmin = QFe[n]->YIndex(PixYmin);
	      jmax = QFe[n]->YIndex(PixYmax);
	      
	      for (i=imin; i<imax; i++)
		{
		  for (j=jmin; j<jmax; j++)
		    {
		      index2 = i + j * QFe[n]->nx;
		      QFe[n]->data[index2] = PixelQFe[m][q];
		    }
		}
	    }
	}
    }
  printf("Finished AdjustQFe, Total Electrons = %.0f\n",TotalElectrons);
  return;
}

double MultiGrid::ElectronQF(int ne)

// Looks up QFe from the QFeLookup table
// With interpolation.
{
  int i = 0;
  double dne, qf, dqf;
  dqf = (QFemax - QFemin) / ((double)NQFe - 1.0);  
  if (ne < 1)
    {
      return 40.0;
    }
  dne = (double)ne;
  if (ne > QFeLookup[0])
    {
      // Ne larger than largest value in lookup table
      // Extrapolate from last two values
      qf = QFemin - dqf / (QFeLookup[0] - QFeLookup[1]) * (dne - QFeLookup[0]);
      //printf("In QFe lookup 0, ne = %d, i = %d, qf = %f\n", ne, i, qf); 
    }
  else
    {
      // Find Ne in the table
      for (i=0; i<NQFe-1; i++)
	{
	  if (ne < QFeLookup[i] && ne > QFeLookup[i+1])
	    {
	      break;
	    }
	}
      if (ne > 1000)  
	{
	  // Ne in table, but larger than 1000
	  // Use linear interpolation
	  qf = QFemin + (double)i * dqf - dqf / (QFeLookup[i] - QFeLookup[i+1]) * (dne - QFeLookup[i]);
	  //printf("In QFe lookup 1, ne = %d, i = %d, qfi = %f, qfi+1 = %f, qf = %f\n", ne, i, QFeLookup[i], QFeLookup[i+1], qf); 
	}
      else
	{
	  // Ne in table, but less than 1000
	  // Use log interpolation
	  // TODO need to add a test for log(0)
	  qf = QFemin + (double)i * dqf - dqf / (log(QFeLookup[i]) - log(QFeLookup[i+1])) * (log(dne) - log(QFeLookup[i]));
	  //printf("In QFe lookup 2, ne = %d, i = %d, qfi = %f, qfi+1 = %f, qf = %f\n", ne, i, QFeLookup[i], QFeLookup[i+1], qf); 
	}
    }
  return max(0.0,min(40.0,qf));
}



void MultiGrid::Setkmins(Array3D** rho, Array3D** phi, Array3D** eps)
{
  int i, j, k, n, index, PixX, Ckmin, Vkmin, NCS, NTaper, NGrid, LastCS = 0, LastC = 0;
  double PixXmin, TaperRatio;
  for (n=0; n<nsteps+1; n++)
    {
      rho[n]->ChannelStopVkmin = 1;
      rho[n]->ChannelStopCkmin = max(rho[n]->ChannelStopVkmin+1, rho[n]->ZIndex(FieldOxide + rho[n]->zmz[rho[n]->ChannelStopVkmin]));
      rho[n]->ChannelVkmin = max(1,rho[n]->ZIndex(FieldOxide / 2.0));
      rho[n]->ChannelCkmin = max(rho[n]->ChannelVkmin+1, rho[n]->ZIndex(GateOxide + rho[n]->zmz[rho[n]->ChannelVkmin]));
      phi[n]->ChannelStopVkmin = rho[n]->ChannelStopVkmin;
      phi[n]->ChannelStopCkmin = rho[n]->ChannelStopCkmin;
      phi[n]->ChannelVkmin = rho[n]->ChannelVkmin;
      phi[n]->ChannelCkmin = rho[n]->ChannelCkmin;  
      //printf("n = %d, CSVkmin = %d, CSCkmin = %d CVkmin = %d, CCkmin = %d\n",n,rho[n]->ChannelStopVkmin,rho[n]->ChannelStopCkmin,rho[n]->ChannelVkmin,rho[n]->ChannelCkmin);
      //fflush(stdout);

      // Count the width in grid of the Channel Stop and taper regions
      NCS = 0; NTaper = 0; NGrid = 0;
      for (i=0; i<phi[n]->nx; i++)
	{
	  PixX = (int)floor((rho[n]->x[i] - PixelRegionLowerLeft[0][0]) / PixelSize);
	  PixXmin = PixelRegionLowerLeft[0][0] + (double)PixX * PixelSize;
	  //printf("n = %d, i = %d, PixX = %d \n",n,i,PixX);
	  //fflush(stdout);

	  if (PixX !=1) continue;
	  if (rho[n]->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho[n]->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0)
	    {
	      // Channel Stop region
	      NCS += 1; NGrid += 1;
	    }
	  else if (rho[n]->x[i] <= PixXmin + ChannelStopWidth/2.0 + FieldOxideTaper)
	    {
	      // Left side taper		      
	      NTaper += 1; NGrid += 1;
	    }
	  else if (rho[n]->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0 - FieldOxideTaper)
	    {
	      // Right side taper		      		      
	      NTaper += 1; NGrid += 1;
	    }
	  else
	    {
	      // Channel region
	      NGrid += 1;
	    }
	}
      //printf("n = %d, NCS = %d, NTaper = %d, NGrid = %d \n",n,NCS,NTaper,NGrid);
      //fflush(stdout);

      for (i=0; i<phi[n]->nx; i++)
	{
	  PixX = (int)floor((rho[n]->x[i] - PixelRegionLowerLeft[0][0]) / PixelSize);
	  PixXmin = PixelRegionLowerLeft[0][0] + (double)PixX * PixelSize;
	  if (rho[n]->x[i] <= PixXmin + ChannelStopWidth/2.0 || rho[n]->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0)
	    {
	      // Channel Stop region
	      TaperRatio = 1.0;
	      LastCS = i;
	    }
	  else if (rho[n]->x[i] <= PixXmin + ChannelStopWidth/2.0 + FieldOxideTaper)
	    {
	      // Left side taper		      
	      TaperRatio = 1.0 - (double)(i - LastCS) / (double)(NTaper/2 + 1);
	    }
	  else if (rho[n]->x[i] >= PixXmin + PixelSize - ChannelStopWidth/2.0 - FieldOxideTaper)
	    {
	      // Right side taper		      		      
	      TaperRatio = (double)(i - LastC) / (double)(NTaper/2 + 1);
	    }
	  else
	    {
	      // Channel region
	      TaperRatio = 0.0;
	      LastC = i;
	    }
	  //if (rho[n]->x[i]>45.0 && rho[n]->x[i]<65.0) printf("n = %d, i = %d, x=%f, TaperRatio = %f\n",n,i,rho[n]->x[i],TaperRatio);
	  //fflush(stdout);
	  Ckmin = rho[n]->ChannelCkmin + (int)round((TaperRatio * (double)(rho[n]->ChannelStopCkmin - rho[n]->ChannelCkmin)));
						    Vkmin = rho[n]->ChannelVkmin + (int)round((TaperRatio * (double)(rho[n]->ChannelStopVkmin - rho[n]->ChannelVkmin)));	  
	  phi[n]->Ckmin[i] = Ckmin;
	  phi[n]->Vkmin[i] = Vkmin;	  
	  rho[n]->Ckmin[i] = Ckmin;
	  rho[n]->Vkmin[i] = Vkmin;
	  rho[n]->TR[i] = TaperRatio;	  	  
	  // Now fill the epsilon arrays
	  for (j=0; j<eps[n]->ny; j++)
	    {
	      for (k=0; k<eps[n]->nz; k++)
		{
		  index = i + j * eps[n]->nx + k * eps[n]->nx * eps[n]->ny;
		  if (k < phi[n]->Ckmin[i]) eps[n]->data[index] = EPSILON_OX / EPSILON_SI;
		  else  eps[n]->data[index] = 1.0;
		}
	    }
	}
    }
  return;
}

void MultiGrid::CountCharges(Array3D** rho, Array3D** elec, Array3D** hole)
{
  // Count total Charge
  int i, j, k, n, index;
  double CellVolume;
  double TotalElectrons = 0.0, TotalHoles = 0.0, TotalCharge = 0.0, ElectronCharge = 0.0, HoleCharge = 0.0;
  for (n=0; n<nsteps+1; n++)
    {
      double RhoChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / (rho[n]->dx * rho[n]->dy);
      TotalElectrons = 0.0; TotalHoles = 0.0; TotalCharge = 0.0; ElectronCharge = 0.0; HoleCharge = 0.0;
      for (i=0; i<rho[n]->nx; i++)
	{
	  for (j=0; j<rho[n]->ny; j++)
	    {
	      for (k=0; k<rho[n]->nz; k++)
		{
		  CellVolume = rho[n]->dx * rho[n]->dy * rho[n]->zw[k];
		  index = i + j * rho[n]->nx + k * rho[n]->nx * rho[n]->ny;
		  TotalCharge += rho[n]->data[index] * CellVolume;
		  if  (k < elec[n]->nz)
		    {
		      TotalHoles += hole[n]->data[index];
		      TotalElectrons += elec[n]->data[index];		      
		      HoleCharge += hole[n]->data[index] * RhoChargeFactor / rho[n]->zw[k] * CellVolume;
		      ElectronCharge += -elec[n]->data[index] * RhoChargeFactor / rho[n]->zw[k] * CellVolume;		      
		    }
		}
	    }
	}
      printf("n = %d, TotalCharge = %.2f, HoleCharge = %.2f, ElectronCharge = %.2f, TotalElectrons = %.1f, TotalHoles = %.1f\n",n,TotalCharge, HoleCharge, ElectronCharge, TotalElectrons, TotalHoles);
    }
  return;
}

void MultiGrid::FillRho(Array3D* rho, Array3D* elec)
{
  //Fill rho with data from electron array
  int i, j, k, index;
  double RhoChargeFactor =  (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / (rho->dx * rho->dy);
  double TotalElectrons = 0.0;
  for (i=0; i<rho->nx; i++)
    { 
      for (j=0; j<rho->ny; j++)
	{
	  for (k=0; k<elec->nz; k++)
	    {
	      index = i + j * rho->nx + k * rho->nx * rho->ny;
	      rho->data[index] += - elec->data[index] * RhoChargeFactor / rho->zw[k];
	      if (elec->data[index] < -1.0E-6)
		{
		  printf("Negative electron count! Something failed!\n");
		}
	      TotalElectrons += elec->data[index];
	    }
	}
    }
  if(VerboseLevel > 1)  printf("Electron charges added into rho. Total electrons=%.1f \n", TotalElectrons);
  return;
}
