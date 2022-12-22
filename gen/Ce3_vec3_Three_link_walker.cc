/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:43 GMT-05:00
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
  double t441;
  double t461;
  double t469;
  double t463;
  double t483;
  double t485;
  double t477;
  double t478;
  double t479;
  double t480;
  double t481;
  double t488;
  double t489;
  double t490;
  double t491;
  double t492;
  double t462;
  double t475;
  double t476;
  double t497;
  double t498;
  double t499;
  double t484;
  double t486;
  double t487;
  double t501;
  double t502;
  double t503;
  t441 = Sin(var1[2]);
  t461 = Cos(var1[3]);
  t469 = Sin(var1[3]);
  t463 = Cos(var1[2]);
  t483 = Cos(var1[4]);
  t485 = Sin(var1[4]);
  t477 = Power(t461,2);
  t478 = -0.5*t477;
  t479 = Power(t469,2);
  t480 = -0.5*t479;
  t481 = t478 + t480;
  t488 = Power(t483,2);
  t489 = -0.5*t488;
  t490 = Power(t485,2);
  t491 = -0.5*t490;
  t492 = t489 + t491;
  t462 = -1.*t461*t441;
  t475 = -1.*t463*t469;
  t476 = t462 + t475;
  t497 = -1.*t463*t461;
  t498 = t441*t469;
  t499 = t497 + t498;
  t484 = -1.*t483*t441;
  t486 = -1.*t463*t485;
  t487 = t484 + t486;
  t501 = -1.*t463*t483;
  t502 = t441*t485;
  t503 = t501 + t502;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(-0.5*(10.*t441 - 5.*t476*t481 - 5.*t487*t492)*var2[0] - 0.5*(10.*t463 - 5.*t481*t499 - 5.*t492*t503)*var2[1])*var2[2];
  p_output1[3]=(2.5*t476*t481*var2[0] + 2.5*t481*t499*var2[1])*var2[2];
  p_output1[4]=(2.5*t487*t492*var2[0] + 2.5*t492*t503*var2[1])*var2[2];
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

#include "Ce3_vec3_Three_link_walker.hh"

namespace SymFunction
{

void Ce3_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
