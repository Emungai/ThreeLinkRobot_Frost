/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:45 GMT-04:00
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
  double t7922;
  double t7923;
  double t7924;
  double t7944;
  double t7964;
  double t7967;
  double t7968;
  double t8000;
  double t8004;
  double t8005;
  double t8006;
  double t8007;
  double t8008;
  double t8012;
  double t8041;
  double t8042;
  double t8043;
  double t8044;
  double t8045;
  double t8046;
  double t8047;
  t7922 = Cos(var1[2]);
  t7923 = Cos(var1[3]);
  t7924 = -1.*t7922*t7923;
  t7944 = Sin(var1[2]);
  t7964 = Sin(var1[3]);
  t7967 = t7944*t7964;
  t7968 = t7924 + t7967;
  t8000 = 1.25*var2[2]*t7968;
  t8004 = 1.25*var2[3]*t7968;
  t8005 = t8000 + t8004;
  t8006 = var2[3]*t8005;
  t8007 = -1.*t7923*t7944;
  t8008 = -1.*t7922*t7964;
  t8012 = t8007 + t8008;
  t8041 = t7923*t7944;
  t8042 = t7922*t7964;
  t8043 = t8041 + t8042;
  t8044 = 1.25*var2[2]*t8043;
  t8045 = 1.25*var2[3]*t8043;
  t8046 = t8044 + t8045;
  t8047 = var2[3]*t8046;
  p_output1[0]=t8006;
  p_output1[1]=t8006;
  p_output1[2]=1.25*t8012*var2[3];
  p_output1[3]=1.25*t8012*var2[2] + 2.5*t8012*var2[3];
  p_output1[4]=t8047;
  p_output1[5]=t8047;
  p_output1[6]=t8004;
  p_output1[7]=t8000 + 2.5*t7968*var2[3];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 8, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_Ce1_vec4_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce1_vec4_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
