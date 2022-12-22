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
  double t7582;
  double t7603;
  double t7583;
  double t7604;
  double t7624;
  double t7625;
  double t7626;
  double t7602;
  double t7610;
  double t7616;
  double t7618;
  double t7620;
  double t7632;
  double t7633;
  double t7634;
  double t7619;
  double t7621;
  double t7622;
  double t7627;
  double t7628;
  double t7629;
  double t7635;
  double t7636;
  double t7637;
  double t7617;
  double t7663;
  double t7664;
  double t7665;
  double t7666;
  double t7667;
  double t7630;
  double t7631;
  double t7642;
  double t7643;
  double t7644;
  double t7645;
  double t7646;
  double t7647;
  double t7648;
  double t7649;
  double t7650;
  double t7623;
  double t7669;
  double t7670;
  double t7671;
  double t7672;
  double t7673;
  double t7638;
  double t7639;
  double t7651;
  double t7652;
  double t7653;
  double t7654;
  double t7655;
  double t7656;
  double t7657;
  double t7658;
  double t7659;
  t7582 = Cos(var1[2]);
  t7603 = Sin(var1[2]);
  t7583 = Cos(var1[3]);
  t7604 = Sin(var1[3]);
  t7624 = -1.*t7583*t7603;
  t7625 = -1.*t7582*t7604;
  t7626 = t7624 + t7625;
  t7602 = -1.*t7582*t7583;
  t7610 = t7603*t7604;
  t7616 = t7602 + t7610;
  t7618 = Cos(var1[4]);
  t7620 = Sin(var1[4]);
  t7632 = -1.*t7618*t7603;
  t7633 = -1.*t7582*t7620;
  t7634 = t7632 + t7633;
  t7619 = -1.*t7582*t7618;
  t7621 = t7603*t7620;
  t7622 = t7619 + t7621;
  t7627 = t7582*t7583;
  t7628 = -1.*t7603*t7604;
  t7629 = t7627 + t7628;
  t7635 = t7582*t7618;
  t7636 = -1.*t7603*t7620;
  t7637 = t7635 + t7636;
  t7617 = -1.25*var2[3]*t7616;
  t7663 = Power(t7583,2);
  t7664 = -0.5*t7663;
  t7665 = Power(t7604,2);
  t7666 = -0.5*t7665;
  t7667 = t7664 + t7666;
  t7630 = -10.*t7626*t7629;
  t7631 = -10.*t7626*t7616;
  t7642 = Power(t7626,2);
  t7643 = -5.*t7642;
  t7644 = t7583*t7603;
  t7645 = t7582*t7604;
  t7646 = t7644 + t7645;
  t7647 = -5.*t7626*t7646;
  t7648 = Power(t7629,2);
  t7649 = -5.*t7648;
  t7650 = -5.*t7629*t7616;
  t7623 = -1.25*var2[4]*t7622;
  t7669 = Power(t7618,2);
  t7670 = -0.5*t7669;
  t7671 = Power(t7620,2);
  t7672 = -0.5*t7671;
  t7673 = t7670 + t7672;
  t7638 = -10.*t7634*t7637;
  t7639 = -10.*t7634*t7622;
  t7651 = Power(t7634,2);
  t7652 = -5.*t7651;
  t7653 = t7618*t7603;
  t7654 = t7582*t7620;
  t7655 = t7653 + t7654;
  t7656 = -5.*t7634*t7655;
  t7657 = Power(t7637,2);
  t7658 = -5.*t7657;
  t7659 = -5.*t7637*t7622;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(t7617 + t7623 - 0.5*(t7643 + t7647 + t7649 + t7650 + t7652 + t7656 + t7658 + t7659)*var2[0] - 0.5*(t7630 + t7631 + t7638 + t7639)*var2[1] - 0.5*(10.*t7582 - 5.*t7616*t7667 - 5.*t7622*t7673)*var2[2]);
  p_output1[3]=var2[1]*(t7617 - 0.5*(t7643 + t7647 + t7649 + t7650)*var2[0] - 0.5*(t7630 + t7631)*var2[1] + 2.5*t7616*t7667*var2[2]);
  p_output1[4]=var2[1]*(t7623 - 0.5*(t7652 + t7656 + t7658 + t7659)*var2[0] - 0.5*(t7638 + t7639)*var2[1] + 2.5*t7622*t7673*var2[2]);
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

#include "Ce3_vec2_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce3_vec2_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
