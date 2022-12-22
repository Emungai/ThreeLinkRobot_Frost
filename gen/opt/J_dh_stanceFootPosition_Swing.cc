/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:45:00 GMT-04:00
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
  double t8632;
  double t8633;
  double t8634;
  double t8638;
  double t8639;
  double t8640;
  double t8641;
  double t8642;
  double t8643;
  double t8644;
  double t8645;
  double t8646;
  double t8647;
  double t8648;
  double t8649;
  double t8650;
  double t8651;
  double t8652;
  double t8653;
  t8632 = Cos(var1[3]);
  t8633 = Sin(var1[2]);
  t8634 = t8632*t8633;
  t8638 = Cos(var1[2]);
  t8639 = Sin(var1[3]);
  t8640 = t8638*t8639;
  t8641 = t8634 + t8640;
  t8642 = var2[2]*t8641;
  t8643 = var2[3]*t8641;
  t8644 = t8642 + t8643;
  t8645 = -1.*t8638*t8632;
  t8646 = t8633*t8639;
  t8647 = t8645 + t8646;
  t8648 = t8638*t8632;
  t8649 = -1.*t8633*t8639;
  t8650 = t8648 + t8649;
  t8651 = var2[2]*t8650;
  t8652 = var2[3]*t8650;
  t8653 = t8651 + t8652;
  p_output1[0]=t8644;
  p_output1[1]=t8644;
  p_output1[2]=1.;
  p_output1[3]=t8647;
  p_output1[4]=t8647;
  p_output1[5]=t8653;
  p_output1[6]=t8653;
  p_output1[7]=1.;
  p_output1[8]=t8641;
  p_output1[9]=t8641;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 10, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_dh_stanceFootPosition_Swing.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_dh_stanceFootPosition_Swing_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
