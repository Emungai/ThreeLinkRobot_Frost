/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:37 GMT-04:00
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
  double t7438;
  double t7444;
  double t7452;
  double t7460;
  double t7461;
  double t7462;
  double t7482;
  double t7504;
  double t7505;
  double t7506;
  double t7512;
  double t7513;
  double t7514;
  double t7523;
  double t7524;
  double t7525;
  double t7498;
  double t7499;
  double t7500;
  double t7501;
  double t7502;
  double t7511;
  double t7518;
  double t7519;
  double t7520;
  double t7521;
  double t7522;
  double t7526;
  double t7527;
  t7438 = Cos(var1[3]);
  t7444 = Sin(var1[2]);
  t7452 = -1.*t7438*t7444;
  t7460 = Cos(var1[2]);
  t7461 = Sin(var1[3]);
  t7462 = -1.*t7460*t7461;
  t7482 = t7452 + t7462;
  t7504 = t7460*t7438;
  t7505 = -1.*t7444*t7461;
  t7506 = t7504 + t7505;
  t7512 = t7438*t7444;
  t7513 = t7460*t7461;
  t7514 = t7512 + t7513;
  t7523 = -1.*t7460*t7438;
  t7524 = t7444*t7461;
  t7525 = t7523 + t7524;
  t7498 = Power(t7438,2);
  t7499 = -0.5*t7498;
  t7500 = Power(t7461,2);
  t7501 = -0.5*t7500;
  t7502 = t7499 + t7501;
  t7511 = 10.*t7482*t7506;
  t7518 = Power(t7482,2);
  t7519 = 5.*t7518;
  t7520 = 5.*t7482*t7514;
  t7521 = Power(t7506,2);
  t7522 = 5.*t7521;
  t7526 = 5.*t7506*t7525;
  t7527 = t7519 + t7520 + t7522 + t7526;
  p_output1[0]=var2[3]*(-0.5*(t7511 + 10.*t7506*t7514)*var2[0] - 0.5*t7527*var2[1] - 2.5*t7482*t7502*var2[2] + 1.25*t7482*var2[3]);
  p_output1[1]=var2[3]*(-0.5*t7527*var2[0] - 0.5*(t7511 + 10.*t7482*t7525)*var2[1] - 2.5*t7502*t7525*var2[2] + 1.25*t7525*var2[3]);
  p_output1[2]=(-2.5*t7482*t7502*var2[0] - 2.5*t7502*t7525*var2[1])*var2[3];
  p_output1[3]=(1.25*t7482*var2[0] + 1.25*t7525*var2[1])*var2[3];
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

#include "Ce2_vec4_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce2_vec4_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
