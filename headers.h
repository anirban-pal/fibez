#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

//INCLUDE GSL HEADERS
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include <cuba.h>
#include <nlopt.h>


typedef std::stringstream sss;

//1D and 2D INTEGRATION parameters 
#define ABS_ERR 0		
#define REL_ERR 1e-4
#define REL_ERR2 1e-12

//Magnitude of random perturbations provided to initial conifguration
#define RTOL 0.0			

//Cuhre integration parameters, see Cuba-Cuhre documentation
#define NVEC 1				
#define EPSREL 1e-12		
#define EPSABS 0
#define VERBOSE 0
#define LAST 4
#define MINEVAL 0
#define MAXEVAL 50000
#define STATEFILE NULL
#define SPIN NULL
#define KEY 0

#define BTYPES 1	//number of cubic bezier types 
double EA[BTYPES] = {0}; //array for storing Axial elastic constants (EA)
double EI[BTYPES] = {0}; //array for storing Bending stiffness constants (EI)
double rho0[BTYPES] = {0}; //array for storing initial mass density of bezier (rho)

#define NBEZ 2e3	//Bezier output resolution in dump file. Determines number of points per bezier in Output.lammpstrj

#define tao_omega 70.0	//Tao integration coupling parameter

#define SVD_factor 1e-4	//factor to check if an singular value from SVD decomp can be treated as 0
#define VERBOSE1 0		//flag to set if intermediate step output is needed
#define C1_flag 1		//enforce C1 continuity if set to 1
#define C2_flag 1		//enforce C2 continuity if set to 1
#define dt 1e-3			//time step
#define NSTEP 1e3		//number of steps

//names of output files
#define DUMP1 "Output_cps.lammpstrj"
#define DUMP2 "Output.lammpstrj"
#define LOG "log.fibez"

#include "potential.h"
#include "util.h"
#include "cbezier.h"
#include "inter.h"
#include "fibnetwork.h"
#include "tao.h"
