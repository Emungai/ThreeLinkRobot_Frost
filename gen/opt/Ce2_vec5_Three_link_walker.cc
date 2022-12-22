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
  double t7497;
  double t7503;
  double t7515;
  double t7516;
  double t7517;
  double t7528;
  double t7529;
  double t7537;
  double t7538;
  double t7539;
  double t7541;
  double t7542;
  double t7543;
  double t7552;
  double t7553;
  double t7554;
  double t7531;
  double t7532;
  double t7533;
  double t7534;
  double t7535;
  double t7540;
  double t7547;
  double t7548;
  double t7549;
  double t7550;
  double t7551;
  double t7555;
  double t7556;
  t7497 = Cos(var1[4]);
  t7503 = Sin(var1[2]);
  t7515 = -1.*t7497*t7503;
  t7516 = Cos(var1[2]);
  t7517 = Sin(var1[4]);
  t7528 = -1.*t7516*t7517;
  t7529 = t7515 + t7528;
  t7537 = t7516*t7497;
  t7538 = -1.*t7503*t7517;
  t7539 = t7537 + t7538;
  t7541 = t7497*t7503;
  t7542 = t7516*t7517;
  t7543 = t7541 + t7542;
  t7552 = -1.*t7516*t7497;
  t7553 = t7503*t7517;
  t7554 = t7552 + t7553;
  t7531 = Power(t7497,2);
  t7532 = -0.5*t7531;
  t7533 = Power(t7517,2);
  t7534 = -0.5*t7533;
  t7535 = t7532 + t7534;
  t7540 = 10.*t7529*t7539;
  t7547 = Power(t7529,2);
  t7548 = 5.*t7547;
  t7549 = 5.*t7529*t7543;
  t7550 = Power(t7539,2);
  t7551 = 5.*t7550;
  t7555 = 5.*t7539*t7554;
  t7556 = t7548 + t7549 + t7551 + t7555;
  p_output1[0]=var2[4]*(-0.5*(t7540 + 10.*t7539*t7543)*var2[0] - 0.5*t7556*var2[1] - 2.5*t7529*t7535*var2[2] + 1.25*t7529*var2[4]);
  p_output1[1]=var2[4]*(-0.5*t7556*var2[0] - 0.5*(t7540 + 10.*t7529*t7554)*var2[1] - 2.5*t7535*t7554*var2[2] + 1.25*t7554*var2[4]);
  p_output1[2]=(-2.5*t7529*t7535*var2[0] - 2.5*t7535*t7554*var2[1])*var2[4];
  p_output1[3]=0;
  p_output1[4]=(1.25*t7529*var2[0] + 1.25*t7554*var2[1])*var2[4];
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

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce2_vec5_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
