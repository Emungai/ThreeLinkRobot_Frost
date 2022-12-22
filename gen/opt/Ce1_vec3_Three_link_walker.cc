/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:34 GMT-04:00
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
  double t7349;
  double t7373;
  double t7350;
  double t7372;
  double t7386;
  double t7388;
  double t7371;
  double t7374;
  double t7375;
  double t7376;
  double t7379;
  double t7380;
  double t7383;
  double t7384;
  double t7387;
  double t7394;
  double t7400;
  double t7401;
  double t7402;
  double t7403;
  double t7404;
  double t7405;
  double t7414;
  double t7415;
  double t7416;
  double t7418;
  double t7419;
  double t7420;
  t7349 = Cos(var1[3]);
  t7373 = Sin(var1[3]);
  t7350 = Sin(var1[2]);
  t7372 = Cos(var1[2]);
  t7386 = Cos(var1[4]);
  t7388 = Sin(var1[4]);
  t7371 = -1.*t7349*t7350;
  t7374 = -1.*t7372*t7373;
  t7375 = t7371 + t7374;
  t7376 = Power(t7349,2);
  t7379 = -0.5*t7376;
  t7380 = Power(t7373,2);
  t7383 = -0.5*t7380;
  t7384 = t7379 + t7383;
  t7387 = -1.*t7386*t7350;
  t7394 = -1.*t7372*t7388;
  t7400 = t7387 + t7394;
  t7401 = Power(t7386,2);
  t7402 = -0.5*t7401;
  t7403 = Power(t7388,2);
  t7404 = -0.5*t7403;
  t7405 = t7402 + t7404;
  t7414 = -1.*t7372*t7349;
  t7415 = t7350*t7373;
  t7416 = t7414 + t7415;
  t7418 = -1.*t7372*t7386;
  t7419 = t7350*t7388;
  t7420 = t7418 + t7419;
  p_output1[0]=var2[2]*(-0.5*(-10.*t7350 + 5.*t7375*t7384 + 5.*t7400*t7405)*var2[2] - 2.5*t7375*t7384*var2[3] - 2.5*t7400*t7405*var2[4]);
  p_output1[1]=var2[2]*(-0.5*(-10.*t7372 + 5.*t7384*t7416 + 5.*t7405*t7420)*var2[2] - 2.5*t7384*t7416*var2[3] - 2.5*t7405*t7420*var2[4]);
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0;
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

#include "Ce1_vec3_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce1_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
