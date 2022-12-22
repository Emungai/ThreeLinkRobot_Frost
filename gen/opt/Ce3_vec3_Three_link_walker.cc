/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:39 GMT-04:00
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
  double t7640;
  double t7660;
  double t7668;
  double t7662;
  double t7682;
  double t7684;
  double t7676;
  double t7677;
  double t7678;
  double t7679;
  double t7680;
  double t7687;
  double t7688;
  double t7689;
  double t7690;
  double t7691;
  double t7661;
  double t7674;
  double t7675;
  double t7696;
  double t7697;
  double t7698;
  double t7683;
  double t7685;
  double t7686;
  double t7700;
  double t7701;
  double t7702;
  t7640 = Sin(var1[2]);
  t7660 = Cos(var1[3]);
  t7668 = Sin(var1[3]);
  t7662 = Cos(var1[2]);
  t7682 = Cos(var1[4]);
  t7684 = Sin(var1[4]);
  t7676 = Power(t7660,2);
  t7677 = -0.5*t7676;
  t7678 = Power(t7668,2);
  t7679 = -0.5*t7678;
  t7680 = t7677 + t7679;
  t7687 = Power(t7682,2);
  t7688 = -0.5*t7687;
  t7689 = Power(t7684,2);
  t7690 = -0.5*t7689;
  t7691 = t7688 + t7690;
  t7661 = -1.*t7660*t7640;
  t7674 = -1.*t7662*t7668;
  t7675 = t7661 + t7674;
  t7696 = -1.*t7662*t7660;
  t7697 = t7640*t7668;
  t7698 = t7696 + t7697;
  t7683 = -1.*t7682*t7640;
  t7685 = -1.*t7662*t7684;
  t7686 = t7683 + t7685;
  t7700 = -1.*t7662*t7682;
  t7701 = t7640*t7684;
  t7702 = t7700 + t7701;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(-0.5*(10.*t7640 - 5.*t7675*t7680 - 5.*t7686*t7691)*var2[0] - 0.5*(10.*t7662 - 5.*t7680*t7698 - 5.*t7691*t7702)*var2[1])*var2[2];
  p_output1[3]=(2.5*t7675*t7680*var2[0] + 2.5*t7680*t7698*var2[1])*var2[2];
  p_output1[4]=(2.5*t7686*t7691*var2[0] + 2.5*t7691*t7702*var2[1])*var2[2];
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

#include "Ce3_vec3_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce3_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
