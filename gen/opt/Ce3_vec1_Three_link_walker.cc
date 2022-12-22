/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:38 GMT-04:00
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
  double t7536;
  double t7545;
  double t7530;
  double t7544;
  double t7546;
  double t7557;
  double t7558;
  double t7566;
  double t7567;
  double t7568;
  double t7560;
  double t7561;
  double t7562;
  double t7563;
  double t7564;
  double t7574;
  double t7575;
  double t7576;
  double t7570;
  double t7571;
  double t7572;
  double t7578;
  double t7579;
  double t7580;
  double t7559;
  double t7605;
  double t7606;
  double t7607;
  double t7608;
  double t7609;
  double t7569;
  double t7573;
  double t7584;
  double t7585;
  double t7586;
  double t7587;
  double t7588;
  double t7589;
  double t7590;
  double t7591;
  double t7592;
  double t7565;
  double t7611;
  double t7612;
  double t7613;
  double t7614;
  double t7615;
  double t7577;
  double t7581;
  double t7593;
  double t7594;
  double t7595;
  double t7596;
  double t7597;
  double t7598;
  double t7599;
  double t7600;
  double t7601;
  t7536 = Sin(var1[2]);
  t7545 = Cos(var1[2]);
  t7530 = Cos(var1[3]);
  t7544 = -1.*t7530*t7536;
  t7546 = Sin(var1[3]);
  t7557 = -1.*t7545*t7546;
  t7558 = t7544 + t7557;
  t7566 = t7545*t7530;
  t7567 = -1.*t7536*t7546;
  t7568 = t7566 + t7567;
  t7560 = Cos(var1[4]);
  t7561 = -1.*t7560*t7536;
  t7562 = Sin(var1[4]);
  t7563 = -1.*t7545*t7562;
  t7564 = t7561 + t7563;
  t7574 = t7545*t7560;
  t7575 = -1.*t7536*t7562;
  t7576 = t7574 + t7575;
  t7570 = t7530*t7536;
  t7571 = t7545*t7546;
  t7572 = t7570 + t7571;
  t7578 = t7560*t7536;
  t7579 = t7545*t7562;
  t7580 = t7578 + t7579;
  t7559 = -1.25*var2[3]*t7558;
  t7605 = Power(t7530,2);
  t7606 = -0.5*t7605;
  t7607 = Power(t7546,2);
  t7608 = -0.5*t7607;
  t7609 = t7606 + t7608;
  t7569 = -10.*t7558*t7568;
  t7573 = -10.*t7572*t7568;
  t7584 = Power(t7558,2);
  t7585 = -5.*t7584;
  t7586 = -5.*t7558*t7572;
  t7587 = Power(t7568,2);
  t7588 = -5.*t7587;
  t7589 = -1.*t7545*t7530;
  t7590 = t7536*t7546;
  t7591 = t7589 + t7590;
  t7592 = -5.*t7568*t7591;
  t7565 = -1.25*var2[4]*t7564;
  t7611 = Power(t7560,2);
  t7612 = -0.5*t7611;
  t7613 = Power(t7562,2);
  t7614 = -0.5*t7613;
  t7615 = t7612 + t7614;
  t7577 = -10.*t7564*t7576;
  t7581 = -10.*t7580*t7576;
  t7593 = Power(t7564,2);
  t7594 = -5.*t7593;
  t7595 = -5.*t7564*t7580;
  t7596 = Power(t7576,2);
  t7597 = -5.*t7596;
  t7598 = -1.*t7545*t7560;
  t7599 = t7536*t7562;
  t7600 = t7598 + t7599;
  t7601 = -5.*t7576*t7600;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(t7559 + t7565 - 0.5*(t7569 + t7573 + t7577 + t7581)*var2[0] - 0.5*(t7585 + t7586 + t7588 + t7592 + t7594 + t7595 + t7597 + t7601)*var2[1] - 0.5*(10.*t7536 - 5.*t7558*t7609 - 5.*t7564*t7615)*var2[2]);
  p_output1[3]=var2[0]*(t7559 - 0.5*(t7569 + t7573)*var2[0] - 0.5*(t7585 + t7586 + t7588 + t7592)*var2[1] + 2.5*t7558*t7609*var2[2]);
  p_output1[4]=var2[0]*(t7565 - 0.5*(t7577 + t7581)*var2[0] - 0.5*(t7594 + t7595 + t7597 + t7601)*var2[1] + 2.5*t7564*t7615*var2[2]);
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

#include "Ce3_vec1_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce3_vec1_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
