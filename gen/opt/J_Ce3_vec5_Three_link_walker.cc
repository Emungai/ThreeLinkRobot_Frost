/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:55 GMT-04:00
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
  double t8571;
  double t8568;
  double t8569;
  double t8572;
  double t8570;
  double t8573;
  double t8574;
  double t8575;
  double t8576;
  double t8577;
  double t8578;
  double t8579;
  double t8580;
  double t8581;
  double t8582;
  double t8583;
  double t8584;
  double t8585;
  double t8586;
  double t8587;
  double t8588;
  double t8589;
  t8571 = Cos(var1[2]);
  t8568 = Cos(var1[4]);
  t8569 = Sin(var1[2]);
  t8572 = Sin(var1[4]);
  t8570 = t8568*t8569;
  t8573 = t8571*t8572;
  t8574 = t8570 + t8573;
  t8575 = -1.25*var2[1]*t8574;
  t8576 = -1.*t8571*t8568;
  t8577 = t8569*t8572;
  t8578 = t8576 + t8577;
  t8579 = -1.25*var2[0]*t8578;
  t8580 = t8575 + t8579;
  t8581 = var2[4]*t8580;
  t8582 = -1.*t8568*t8569;
  t8583 = -1.*t8571*t8572;
  t8584 = t8582 + t8583;
  t8585 = -1.25*var2[4]*t8584;
  t8586 = -1.25*var2[4]*t8578;
  t8587 = -1.25*var2[0]*t8584;
  t8588 = -1.25*var2[1]*t8578;
  t8589 = t8587 + t8588;
  p_output1[0]=t8581;
  p_output1[1]=t8581;
  p_output1[2]=t8585;
  p_output1[3]=t8586;
  p_output1[4]=t8589;
  p_output1[5]=t8581;
  p_output1[6]=t8581;
  p_output1[7]=t8585;
  p_output1[8]=t8586;
  p_output1[9]=t8589;
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

#include "J_Ce3_vec5_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce3_vec5_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
