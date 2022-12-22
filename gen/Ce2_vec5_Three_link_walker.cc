/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:37 GMT-05:00
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
  double t298;
  double t304;
  double t316;
  double t317;
  double t318;
  double t329;
  double t330;
  double t338;
  double t339;
  double t340;
  double t342;
  double t343;
  double t344;
  double t353;
  double t354;
  double t355;
  double t332;
  double t333;
  double t334;
  double t335;
  double t336;
  double t341;
  double t348;
  double t349;
  double t350;
  double t351;
  double t352;
  double t356;
  double t357;
  t298 = Cos(var1[4]);
  t304 = Sin(var1[2]);
  t316 = -1.*t298*t304;
  t317 = Cos(var1[2]);
  t318 = Sin(var1[4]);
  t329 = -1.*t317*t318;
  t330 = t316 + t329;
  t338 = t317*t298;
  t339 = -1.*t304*t318;
  t340 = t338 + t339;
  t342 = t298*t304;
  t343 = t317*t318;
  t344 = t342 + t343;
  t353 = -1.*t317*t298;
  t354 = t304*t318;
  t355 = t353 + t354;
  t332 = Power(t298,2);
  t333 = -0.5*t332;
  t334 = Power(t318,2);
  t335 = -0.5*t334;
  t336 = t333 + t335;
  t341 = 10.*t330*t340;
  t348 = Power(t330,2);
  t349 = 5.*t348;
  t350 = 5.*t330*t344;
  t351 = Power(t340,2);
  t352 = 5.*t351;
  t356 = 5.*t340*t355;
  t357 = t349 + t350 + t352 + t356;
  p_output1[0]=var2[4]*(-0.5*(t341 + 10.*t340*t344)*var2[0] - 0.5*t357*var2[1] - 2.5*t330*t336*var2[2] + 1.25*t330*var2[4]);
  p_output1[1]=var2[4]*(-0.5*t357*var2[0] - 0.5*(t341 + 10.*t330*t355)*var2[1] - 2.5*t336*t355*var2[2] + 1.25*t355*var2[4]);
  p_output1[2]=(-2.5*t330*t336*var2[0] - 2.5*t336*t355*var2[1])*var2[4];
  p_output1[3]=0;
  p_output1[4]=(1.25*t330*var2[0] + 1.25*t355*var2[1])*var2[4];
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

#include "Ce2_vec5_Three_link_walker.hh"

namespace SymFunction
{

void Ce2_vec5_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
