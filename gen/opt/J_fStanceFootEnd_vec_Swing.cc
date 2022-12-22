/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:33:05 GMT-04:00
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
  double t6787;
  double t6784;
  double t6785;
  double t6791;
  double t6786;
  double t6792;
  double t6793;
  double t6794;
  double t6797;
  double t6798;
  double t6799;
  double t6800;
  double t6801;
  double t6805;
  double t6806;
  double t6807;
  double t6808;
  double t6809;
  double t6810;
  t6787 = Cos(var1[2]);
  t6784 = Cos(var1[3]);
  t6785 = Sin(var1[2]);
  t6791 = Sin(var1[3]);
  t6786 = -1.*t6784*t6785;
  t6792 = -1.*t6787*t6791;
  t6793 = t6786 + t6792;
  t6794 = var2[0]*t6793;
  t6797 = t6787*t6784;
  t6798 = -1.*t6785*t6791;
  t6799 = t6797 + t6798;
  t6800 = var2[2]*t6799;
  t6801 = t6794 + t6800;
  t6805 = var2[2]*t6793;
  t6806 = -1.*t6787*t6784;
  t6807 = t6785*t6791;
  t6808 = t6806 + t6807;
  t6809 = var2[0]*t6808;
  t6810 = t6805 + t6809;
  p_output1[0]=t6801;
  p_output1[1]=t6801;
  p_output1[2]=t6799;
  p_output1[3]=t6784*t6785 + t6787*t6791;
  p_output1[4]=t6810;
  p_output1[5]=t6810;
  p_output1[6]=t6793;
  p_output1[7]=t6799;
  p_output1[8]=-1.*Power(t6784,2) - 1.*Power(t6791,2);
  p_output1[9]=-1.;
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
    ( !(mrows == 3 && ncols == 1) && 
      !(mrows == 1 && ncols == 3))) 
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

#include "J_fStanceFootEnd_vec_Swing.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_fStanceFootEnd_vec_Swing_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
