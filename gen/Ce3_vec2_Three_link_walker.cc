/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:42 GMT-05:00
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
  double t383;
  double t404;
  double t384;
  double t405;
  double t425;
  double t426;
  double t427;
  double t403;
  double t411;
  double t417;
  double t419;
  double t421;
  double t433;
  double t434;
  double t435;
  double t420;
  double t422;
  double t423;
  double t428;
  double t429;
  double t430;
  double t436;
  double t437;
  double t438;
  double t418;
  double t464;
  double t465;
  double t466;
  double t467;
  double t468;
  double t431;
  double t432;
  double t443;
  double t444;
  double t445;
  double t446;
  double t447;
  double t448;
  double t449;
  double t450;
  double t451;
  double t424;
  double t470;
  double t471;
  double t472;
  double t473;
  double t474;
  double t439;
  double t440;
  double t452;
  double t453;
  double t454;
  double t455;
  double t456;
  double t457;
  double t458;
  double t459;
  double t460;
  t383 = Cos(var1[2]);
  t404 = Sin(var1[2]);
  t384 = Cos(var1[3]);
  t405 = Sin(var1[3]);
  t425 = -1.*t384*t404;
  t426 = -1.*t383*t405;
  t427 = t425 + t426;
  t403 = -1.*t383*t384;
  t411 = t404*t405;
  t417 = t403 + t411;
  t419 = Cos(var1[4]);
  t421 = Sin(var1[4]);
  t433 = -1.*t419*t404;
  t434 = -1.*t383*t421;
  t435 = t433 + t434;
  t420 = -1.*t383*t419;
  t422 = t404*t421;
  t423 = t420 + t422;
  t428 = t383*t384;
  t429 = -1.*t404*t405;
  t430 = t428 + t429;
  t436 = t383*t419;
  t437 = -1.*t404*t421;
  t438 = t436 + t437;
  t418 = -1.25*var2[3]*t417;
  t464 = Power(t384,2);
  t465 = -0.5*t464;
  t466 = Power(t405,2);
  t467 = -0.5*t466;
  t468 = t465 + t467;
  t431 = -10.*t427*t430;
  t432 = -10.*t427*t417;
  t443 = Power(t427,2);
  t444 = -5.*t443;
  t445 = t384*t404;
  t446 = t383*t405;
  t447 = t445 + t446;
  t448 = -5.*t427*t447;
  t449 = Power(t430,2);
  t450 = -5.*t449;
  t451 = -5.*t430*t417;
  t424 = -1.25*var2[4]*t423;
  t470 = Power(t419,2);
  t471 = -0.5*t470;
  t472 = Power(t421,2);
  t473 = -0.5*t472;
  t474 = t471 + t473;
  t439 = -10.*t435*t438;
  t440 = -10.*t435*t423;
  t452 = Power(t435,2);
  t453 = -5.*t452;
  t454 = t419*t404;
  t455 = t383*t421;
  t456 = t454 + t455;
  t457 = -5.*t435*t456;
  t458 = Power(t438,2);
  t459 = -5.*t458;
  t460 = -5.*t438*t423;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(t418 + t424 - 0.5*(t444 + t448 + t450 + t451 + t453 + t457 + t459 + t460)*var2[0] - 0.5*(t431 + t432 + t439 + t440)*var2[1] - 0.5*(10.*t383 - 5.*t417*t468 - 5.*t423*t474)*var2[2]);
  p_output1[3]=var2[1]*(t418 - 0.5*(t444 + t448 + t450 + t451)*var2[0] - 0.5*(t431 + t432)*var2[1] + 2.5*t417*t468*var2[2]);
  p_output1[4]=var2[1]*(t424 - 0.5*(t453 + t457 + t459 + t460)*var2[0] - 0.5*(t439 + t440)*var2[1] + 2.5*t423*t474*var2[2]);
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

#include "Ce3_vec2_Three_link_walker.hh"

namespace SymFunction
{

void Ce3_vec2_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
