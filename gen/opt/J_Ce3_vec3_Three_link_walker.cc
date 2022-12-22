/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:53 GMT-04:00
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
  double t3230;
  double t3622;
  double t8440;
  double t8439;
  double t8471;
  double t8504;
  double t8455;
  double t8461;
  double t8467;
  double t8468;
  double t8469;
  double t8511;
  double t8512;
  double t8513;
  double t8517;
  double t8518;
  double t4665;
  double t8453;
  double t8454;
  double t8523;
  double t8525;
  double t8526;
  double t8497;
  double t8505;
  double t8506;
  double t8528;
  double t8529;
  double t8530;
  double t8522;
  double t8527;
  double t8531;
  double t8532;
  double t8544;
  double t8545;
  double t8546;
  double t8547;
  double t8548;
  double t8549;
  double t8550;
  double t8551;
  double t8552;
  double t8553;
  double t8536;
  double t8537;
  double t8538;
  double t8539;
  double t8540;
  double t8541;
  double t8542;
  double t8543;
  t3230 = Sin(var1[2]);
  t3622 = Cos(var1[3]);
  t8440 = Sin(var1[3]);
  t8439 = Cos(var1[2]);
  t8471 = Cos(var1[4]);
  t8504 = Sin(var1[4]);
  t8455 = Power(t3622,2);
  t8461 = -0.5*t8455;
  t8467 = Power(t8440,2);
  t8468 = -0.5*t8467;
  t8469 = t8461 + t8468;
  t8511 = Power(t8471,2);
  t8512 = -0.5*t8511;
  t8513 = Power(t8504,2);
  t8517 = -0.5*t8513;
  t8518 = t8512 + t8517;
  t4665 = t3622*t3230;
  t8453 = t8439*t8440;
  t8454 = t4665 + t8453;
  t8523 = -1.*t8439*t3622;
  t8525 = t3230*t8440;
  t8526 = t8523 + t8525;
  t8497 = t8471*t3230;
  t8505 = t8439*t8504;
  t8506 = t8497 + t8505;
  t8528 = -1.*t8439*t8471;
  t8529 = t3230*t8504;
  t8530 = t8528 + t8529;
  t8522 = 10.*t8439;
  t8527 = -5.*t8526*t8469;
  t8531 = -5.*t8530*t8518;
  t8532 = t8522 + t8527 + t8531;
  t8544 = 10.*t3230;
  t8545 = -1.*t3622*t3230;
  t8546 = -1.*t8439*t8440;
  t8547 = t8545 + t8546;
  t8548 = -5.*t8547*t8469;
  t8549 = -1.*t8471*t3230;
  t8550 = -1.*t8439*t8504;
  t8551 = t8549 + t8550;
  t8552 = -5.*t8551*t8518;
  t8553 = t8544 + t8548 + t8552;
  t8536 = 2.5*var2[1]*t8454*t8469;
  t8537 = 2.5*var2[0]*t8526*t8469;
  t8538 = t8536 + t8537;
  t8539 = var2[2]*t8538;
  t8540 = 2.5*var2[1]*t8506*t8518;
  t8541 = 2.5*var2[0]*t8530*t8518;
  t8542 = t8540 + t8541;
  t8543 = var2[2]*t8542;
  p_output1[0]=(-0.5*t8532*var2[0] - 0.5*(-10.*t3230 - 5.*t8454*t8469 - 5.*t8506*t8518)*var2[1])*var2[2];
  p_output1[1]=t8539;
  p_output1[2]=t8543;
  p_output1[3]=-0.5*t8553*var2[2];
  p_output1[4]=-0.5*t8532*var2[2];
  p_output1[5]=-0.5*t8553*var2[0] - 0.5*t8532*var2[1];
  p_output1[6]=t8539;
  p_output1[7]=t8539;
  p_output1[8]=2.5*t8469*t8547*var2[2];
  p_output1[9]=2.5*t8469*t8526*var2[2];
  p_output1[10]=2.5*t8469*t8547*var2[0] + 2.5*t8469*t8526*var2[1];
  p_output1[11]=t8543;
  p_output1[12]=t8543;
  p_output1[13]=2.5*t8518*t8551*var2[2];
  p_output1[14]=2.5*t8518*t8530*var2[2];
  p_output1[15]=2.5*t8518*t8551*var2[0] + 2.5*t8518*t8530*var2[1];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 16, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_Ce3_vec3_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce3_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
