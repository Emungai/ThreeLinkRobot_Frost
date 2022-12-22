/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:46 GMT-04:00
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
  double t3041;
  double t3042;
  double t3047;
  double t8016;
  double t8018;
  double t8019;
  double t8040;
  double t8048;
  double t8049;
  double t8054;
  double t8055;
  double t8056;
  double t8057;
  double t8058;
  double t8072;
  double t8073;
  double t8075;
  double t8076;
  double t8077;
  double t8078;
  double t8079;
  t3041 = Cos(var1[2]);
  t3042 = Cos(var1[4]);
  t3047 = -1.*t3041*t3042;
  t8016 = Sin(var1[2]);
  t8018 = Sin(var1[4]);
  t8019 = t8016*t8018;
  t8040 = t3047 + t8019;
  t8048 = 1.25*var2[2]*t8040;
  t8049 = 1.25*var2[4]*t8040;
  t8054 = t8048 + t8049;
  t8055 = var2[4]*t8054;
  t8056 = -1.*t3042*t8016;
  t8057 = -1.*t3041*t8018;
  t8058 = t8056 + t8057;
  t8072 = t3042*t8016;
  t8073 = t3041*t8018;
  t8075 = t8072 + t8073;
  t8076 = 1.25*var2[2]*t8075;
  t8077 = 1.25*var2[4]*t8075;
  t8078 = t8076 + t8077;
  t8079 = var2[4]*t8078;
  p_output1[0]=t8055;
  p_output1[1]=t8055;
  p_output1[2]=1.25*t8058*var2[4];
  p_output1[3]=1.25*t8058*var2[2] + 2.5*t8058*var2[4];
  p_output1[4]=t8079;
  p_output1[5]=t8079;
  p_output1[6]=t8049;
  p_output1[7]=t8048 + 2.5*t8040*var2[4];
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

#include "J_Ce1_vec5_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce1_vec5_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
