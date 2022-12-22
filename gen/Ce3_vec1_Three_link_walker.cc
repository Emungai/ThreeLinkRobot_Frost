/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:39 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
/**
 * Copied from Wolfram Mathematica C Definitions file mdefs.hpp
 * Changed marcos to inline functions (Eric Cousineau)
 */
inline double Power(double x, double y) { return pow(x, y); }
inline double Sqrt(double x) { return sqrt(x); }

inline double Abs(double x) { return fabs(x); }

inline double Exp(double x) { return exp(x); }
inline double Log(double x) { return log(x); }

inline double Sin(double x) { return sin(x); }
inline double Cos(double x) { return cos(x); }
inline double Tan(double x) { return tan(x); }

inline double ArcSin(double x) { return asin(x); }
inline double ArcCos(double x) { return acos(x); }
inline double ArcTan(double x) { return atan(x); }

/* update ArcTan function to use atan2 instead. */
inline double ArcTan(double x, double y) { return atan2(y,x); }

inline double Sinh(double x) { return sinh(x); }
inline double Cosh(double x) { return cosh(x); }
inline double Tanh(double x) { return tanh(x); }

const double E	= 2.71828182845904523536029;
const double Pi = 3.14159265358979323846264;
const double Degree = 0.01745329251994329576924;


#endif

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t337;
  double t346;
  double t331;
  double t345;
  double t347;
  double t358;
  double t359;
  double t367;
  double t368;
  double t369;
  double t361;
  double t362;
  double t363;
  double t364;
  double t365;
  double t375;
  double t376;
  double t377;
  double t371;
  double t372;
  double t373;
  double t379;
  double t380;
  double t381;
  double t360;
  double t406;
  double t407;
  double t408;
  double t409;
  double t410;
  double t370;
  double t374;
  double t385;
  double t386;
  double t387;
  double t388;
  double t389;
  double t390;
  double t391;
  double t392;
  double t393;
  double t366;
  double t412;
  double t413;
  double t414;
  double t415;
  double t416;
  double t378;
  double t382;
  double t394;
  double t395;
  double t396;
  double t397;
  double t398;
  double t399;
  double t400;
  double t401;
  double t402;
  t337 = Sin(var1[2]);
  t346 = Cos(var1[2]);
  t331 = Cos(var1[3]);
  t345 = -1.*t331*t337;
  t347 = Sin(var1[3]);
  t358 = -1.*t346*t347;
  t359 = t345 + t358;
  t367 = t346*t331;
  t368 = -1.*t337*t347;
  t369 = t367 + t368;
  t361 = Cos(var1[4]);
  t362 = -1.*t361*t337;
  t363 = Sin(var1[4]);
  t364 = -1.*t346*t363;
  t365 = t362 + t364;
  t375 = t346*t361;
  t376 = -1.*t337*t363;
  t377 = t375 + t376;
  t371 = t331*t337;
  t372 = t346*t347;
  t373 = t371 + t372;
  t379 = t361*t337;
  t380 = t346*t363;
  t381 = t379 + t380;
  t360 = -1.25*var2[3]*t359;
  t406 = Power(t331,2);
  t407 = -0.5*t406;
  t408 = Power(t347,2);
  t409 = -0.5*t408;
  t410 = t407 + t409;
  t370 = -10.*t359*t369;
  t374 = -10.*t373*t369;
  t385 = Power(t359,2);
  t386 = -5.*t385;
  t387 = -5.*t359*t373;
  t388 = Power(t369,2);
  t389 = -5.*t388;
  t390 = -1.*t346*t331;
  t391 = t337*t347;
  t392 = t390 + t391;
  t393 = -5.*t369*t392;
  t366 = -1.25*var2[4]*t365;
  t412 = Power(t361,2);
  t413 = -0.5*t412;
  t414 = Power(t363,2);
  t415 = -0.5*t414;
  t416 = t413 + t415;
  t378 = -10.*t365*t377;
  t382 = -10.*t381*t377;
  t394 = Power(t365,2);
  t395 = -5.*t394;
  t396 = -5.*t365*t381;
  t397 = Power(t377,2);
  t398 = -5.*t397;
  t399 = -1.*t346*t361;
  t400 = t337*t363;
  t401 = t399 + t400;
  t402 = -5.*t377*t401;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(t360 + t366 - 0.5*(t370 + t374 + t378 + t382)*var2[0] - 0.5*(t386 + t387 + t389 + t393 + t395 + t396 + t398 + t402)*var2[1] - 0.5*(10.*t337 - 5.*t359*t410 - 5.*t365*t416)*var2[2]);
  p_output1[3]=var2[0]*(t360 - 0.5*(t370 + t374)*var2[0] - 0.5*(t386 + t387 + t389 + t393)*var2[1] + 2.5*t359*t410*var2[2]);
  p_output1[4]=var2[0]*(t366 - 0.5*(t378 + t382)*var2[0] - 0.5*(t395 + t396 + t398 + t402)*var2[1] + 2.5*t365*t416*var2[2]);
}



#ifdef MATLAB_MEX_FILE

#include "mex.h"
/*
 * Main function
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  size_t mrows, ncols;

  double *var1,*var2;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 2)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Two input(s) required (var1,var2).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 5, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "Ce3_vec1_Three_link_walker.hh"

namespace SymFunction
{

void Ce3_vec1_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
