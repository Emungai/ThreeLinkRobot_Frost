/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:54 GMT-04:00
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
  double t8520;
  double t3616;
  double t8470;
  double t8521;
  double t8519;
  double t8533;
  double t8534;
  double t8535;
  double t8554;
  double t8555;
  double t8556;
  double t8557;
  double t8558;
  double t8559;
  double t8560;
  double t8561;
  double t8562;
  double t8563;
  double t8564;
  double t8565;
  double t8566;
  double t8567;
  t8520 = Cos(var1[2]);
  t3616 = Cos(var1[3]);
  t8470 = Sin(var1[2]);
  t8521 = Sin(var1[3]);
  t8519 = t3616*t8470;
  t8533 = t8520*t8521;
  t8534 = t8519 + t8533;
  t8535 = -1.25*var2[1]*t8534;
  t8554 = -1.*t8520*t3616;
  t8555 = t8470*t8521;
  t8556 = t8554 + t8555;
  t8557 = -1.25*var2[0]*t8556;
  t8558 = t8535 + t8557;
  t8559 = var2[3]*t8558;
  t8560 = -1.*t3616*t8470;
  t8561 = -1.*t8520*t8521;
  t8562 = t8560 + t8561;
  t8563 = -1.25*var2[3]*t8562;
  t8564 = -1.25*var2[3]*t8556;
  t8565 = -1.25*var2[0]*t8562;
  t8566 = -1.25*var2[1]*t8556;
  t8567 = t8565 + t8566;
  p_output1[0]=t8559;
  p_output1[1]=t8559;
  p_output1[2]=t8563;
  p_output1[3]=t8564;
  p_output1[4]=t8567;
  p_output1[5]=t8559;
  p_output1[6]=t8559;
  p_output1[7]=t8563;
  p_output1[8]=t8564;
  p_output1[9]=t8567;
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

#include "J_Ce3_vec4_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce3_vec4_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
