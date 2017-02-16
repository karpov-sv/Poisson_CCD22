/*
  ------------------------------------------------------------------------------
  Author: Craig Lage, UC Davis
  Date: Sep 30, 2016

  Standalone cpp Poisson solver

*/

//****************** multigrid.h **************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/iter_find.hpp>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace std;

#include <globals.h>
#include <fileio.h>
#include <hdf5write.h>
#include <hdf5read.h>
#include <array3d.h>
#include <array2d.h>
#include <polygon.h>

class MultiGrid
{
 public:

  double w;		// Successive Over-Relaxation factor
  int ncycle;		// Number of SOR cycles at each resolution
  int ncycle_loop;	// Number of SOR cycles in repeated loops  
  int iterations;	// Number of VCycles
  double NZExp;		// Non-linear z axis exponent
  
  int ScaleFactor;       // Power of 2 that sets the grid size
  // ScaleFactor = 1 means grid size is 1.0 micron, 96 grids in the z-direction
  double PixelSize;      // Pixel size in microns
  double SensorThickness; // Thickness of the sensor in microns
  int GridsPerPixel;     // Number of grids per pixel at ScaleFactor = 1
  double GridSpacing;    // Grid size in microns
  int Nx;                // Number of grids in x at ScaleFactor = 1
  int Ny;                // Number of grids in y at ScaleFactor = 1
  int Nz;                // Number of grids in z at ScaleFactor = 1
  int Nzelec;            // Number of grids in z in electron arrays at ScaleFactor = 1
  double Xmin;
  double Xmax;
  double Ymin;
  double Ymax;
  double Zmin;
  double Zmax;
  double* SimulationRegionLowerLeft;
  int nsteps;            // Number of steps in Multigrid
  int XBCType;            // Free or periodic BC in X-direction
  int YBCType;            // Free or periodic BC in Y-direction

  // Voltages and Charges
  double Ni;            // Intrinsic carrier concentration at operating temperature
  double ktq;           // kT/q
  double Vbb;		// Back bias
  double Vparallel_lo;	// Parallel Low Voltage
  double Vparallel_hi;	// Parallel High Voltage
  double Vserial_lo;	// Serial Low Voltage
  double Vaverage;       // Average voltage on bottom
  double Vchannelstop;   //  Voltage of undepleted channel stop
  double Vscupper;       // Scupper voltage

  int Channelkmin;              // Bottom of channel region doping
  int Channelkmax;              // Top of channel region doping
  int Troughkmin;              // Bottom of trough region doping
  int Troughkmax;              // Top of trough region doping
  int ChannelStopkmin;          // Bottom of channel stop region doping
  int ChannelStopkmax;          // Top of channel stop region doping
  int ChannelProfile;           // 0 = Square well, 1 = Gaussian
  int ChannelStopProfile;       // 0 = Square well, 1 = Gaussian
  double BackgroundDoping; 	// Background doping
  double ChannelStopDoping;	// Channel Stop doping
  double ChannelStopDepth;     	// Channel stop depth in microns
  double ChannelStopWidth;     	// Channel stop width in microns
  double ChannelStopCharge;
  double ChannelStopPeak;       // Depth of peak of channel stop implant below silicon surface in microns
  double CSChargeDepth;
  double ChannelDoping;		// Channel doping
  double ChannelDepth;		// Channel depth in microns
  double ChannelPeak;           // Depth of peak of channel implant below silicon surface in microns
  int TroughImplant;            // 1 = Uses Trough implant; 0 = No Trough implant
  double TroughDoping;          // Trough doping in cm^-2
  int TroughProfile;            // 0 = Square profile, 1 = Gaussian profile
  double TroughDepth;           // Trough depth in microns
  double TroughWidth;           // Trough width in microns
  double TroughPeak;           // Depth of peak of trough implant below silicon surface in microns
  double GateOxide;             // Gate oxide thickness in microns
  double Gox_effective;         // Effective gate oxide thickness in microns (includes effect of discretization)
  double FieldOxide;            // Field oxide thickness in microns
  double FieldOxideTaper;       // Field oxide taper width in microns  
  double Delta_Z;               // Z-axis shift due to difference in dielectric constant between Si and SiO2
  int UndepletedChannelStop;	// 0 = No undepleted Region, 1 = undepleted Region
  double HoleConvergenceVoltage;  // Voltage to which we converge the hole calculation

// Pixel Regions

  int NumberofPixelRegions;	  	  // 
  double** PixelRegionLowerLeft;	  //
  double** PixelRegionUpperRight;	  //
  int* NumberofFilledWells;		  //
  int** CollectedCharge;		  // Collected charge in e-
  double*** FilledPixelCoords;            // (x,y) coords of pixel center
  double** ElectronCount;                 // Number of electrons in each filled pixel
  double** PixelQFe;                      // QFe in each filled pixel

  // QFe Look-up table
  int BuildQFeLookup;        // 0 - No QFe lookup; 1 - Build QFeLookup
  int NQFe;                  // Number of entries in QFELookup table             
  double QFemin;             // Minimum QFe in table
  double QFemax;             // Maximum QFe in table
  double* QFeLookup;         // QFeLookup table
  
  int CollectingPhases;                   // Number of collecting phases

  // Added dipoles

  int NumberofDipoles;                  // Number of dipoles
  double** DipoleCoords;                // X,Y coords of dipole center
  double* DipoleZLocation;               // Z location of dipole charge
  int* DipoleCharge;                     // Number of electrons in dipole
  int NumberofDipoleImages;             // Number of image pairs to calculate

// Constant Voltage Regions

  int NumberofFixedRegions;
  double** FixedRegionLowerLeft;
  double** FixedRegionUpperRight;
  double* FixedRegionVoltage;
  double* FixedRegionDoping;
  double* FixedRegionChargeDepth;
  int* FixedRegionBCType;

  // Pixel Boundary Tests

  int PixelBoundaryTestType;
  int LogEField;
  int LogPixelPaths;
  int PixelAreas;
  double* PixelBoundaryLowerLeft;
  double* PixelBoundaryUpperRight;
  double* PixelBoundaryStepSize;

  int PixelBoundaryNx;
  int PixelBoundaryNy;
  int NumVertices;
  double ElectronZ0Area;
  double ElectronZ0Fill;
  double CCDTemperature;
  double DiffMultiplier;
  int NumDiffSteps;
  int EquilibrateSteps;
  int BottomSteps;
  int ElectronAccumulation;      
  int SaturationModel;
  int NumElec;
  int NumSteps;
  double Sigmax;
  double Sigmay;
  double Xoffset;
  double Yoffset;

  string outputfilebase; // Output filename base
  string outputfiledir; // Output filename directory
  string PhotonList; // Photon list filename
  int NumPhotons;  //Number of photons in list
  double* PhotonListx;
  double* PhotonListy;
  double* PhotonListdxdz;
  double* PhotonListdydz;
  double* PhotonListlambda;  

  int VerboseLevel;
  int SaveData;
  int SaveElec;

// Continuation info

  int Continuation;
  int LastContinuationStep;  


  string FilterBand;    // One of "u", "g", "r", "i", "z", "y".
  int FilterIndex;      // 0=u, 1=g, 2=r, 3=i, 4=z, 5=y
  string FilterFile;    // location of tabulated per-band depth probabilities
  static const int n_band = 6, n_filter_cdf = 5000;
  double filter_cdf[n_band * n_filter_cdf];

  Array3D** phi;      // Phi arrays
  Array3D** rho;      // Rho arrays
  Array3D** E;        // Electric field
  Array3D** elec;      // Number of stored electrons
  Array3D** hole;      // Number of mobile holes
  Array3D** eps;      // Dielectric constant array    
  Array2D** BCType;      // BCType - 0->fixed potential; 1->Enormal = 0
  Array2D** QFe;         // Electron Quasi-Fermi level
  Array2D** QFh;         // Hole Quasi-Fermi level  
  double qfe;
  double qfh;
  MultiGrid() {};
  MultiGrid(string);
  ~MultiGrid();

  void ReadConfigurationFile(string);
  void ReadPhotonList(string, string);
  void BuildArrays(Array3D**, Array3D**, Array3D**, Array3D**, Array3D**, Array3D**, Array2D**, Array2D**, Array2D**);
  void SaveGrid();
  void SetInitialVoltages(Array3D*, Array2D*);
  void SetFixedCharges(Array3D*, Array2D*);
  double SOR_Inner(Array3D*, Array3D*, Array3D*, Array3D*, Array3D*, Array2D*, Array2D*, Array2D*, double);
  double Error_Inner(Array3D*, Array3D*, Array3D*, Array3D*, Array3D*);
  void Restrict(Array3D*, Array3D*, Array3D*, Array3D*, Array2D*, Array2D*);
  void Prolongate(Array3D*, Array3D*, Array2D*);
  void VCycle_Inner(Array3D**, Array3D**, Array3D**, Array3D**, Array3D**, Array2D**, Array2D**, Array2D**, double, int, int);
  void VCycle_Zero(Array3D**, Array3D**, Array3D**, Array3D**, Array3D**, Array2D**, Array2D**, Array2D**, double, int);  
  void WriteOutputFile(string, string, string, Array3D*);
  void ReadOutputFile(string, string, string, Array3D*);
  void Gradient(Array3D*, Array3D**);
  double GetElectronInitialZ();
  void Trace(double*, int, bool, double, ofstream&);
  void TraceSpot(int);
  void TraceList(int);      
  void TraceGrid(int);
  void TraceRegion(int);
  void FindEdge(double*, double, ofstream&);
  void FindCorner(double*, double*, ofstream&);
  void CalculatePixelAreas(int);
  double mu_Si (double, double);
  void Set_QFh(Array2D**, double);
  void Adjust_QFe(Array2D**, Array3D*, Array3D*);
  double ElectronQF(int);
  void WriteQFeLookup(string, string, string);
  void WriteCollectedCharge(string, string, string);  
  void ReadQFeLookup(string, string, string);
  void Setkmins(Array3D**, Array3D**, Array3D**);
  void CountCharges(Array3D**, Array3D**, Array3D**);
  void FillRho(Array3D*, Array3D*);
};
