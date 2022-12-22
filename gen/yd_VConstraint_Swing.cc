/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 16:38:25 GMT-04:00
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
static void output1(double *p_output1,const double *var1,const double *var2,const double *var3)
{
  double t1061;
  double t1062;
  double t1064;
  double t1069;
  double t1070;
  double t1071;
  double t1063;
  double t1065;
  double t1067;
  double t1068;
  double t1073;
  double t1074;
  double t1075;
  double t1077;
  double t1078;
  double t1079;
  double t1081;
  double t1083;
  t1061 = -1.*var3[0];
  t1062 = t1061 + var3[1];
  t1064 = t1061 + var1[2] + var1[3];
  t1069 = 1/t1062;
  t1070 = -1.*t1069*t1064;
  t1071 = 1. + t1070;
  t1063 = Power(t1062,-5);
  t1065 = Power(t1064,5);
  t1067 = Power(t1062,-4);
  t1068 = Power(t1064,4);
  t1073 = Power(t1062,-3);
  t1074 = Power(t1064,3);
  t1075 = Power(t1071,2);
  t1077 = Power(t1062,-2);
  t1078 = Power(t1064,2);
  t1079 = Power(t1071,3);
  t1081 = Power(t1071,4);
  t1083 = Power(t1071,5);
  p_output1[0]=t1083*var2[0] + 5.*t1064*t1069*t1081*var2[2] + 10.*t1077*t1078*t1079*var2[4] + 10.*t1073*t1074*t1075*var2[6] + 5.*t1067*t1068*t1071*var2[8] + t1063*t1065*var2[10];
  p_output1[1]=t1083*var2[1] + 5.*t1064*t1069*t1081*var2[3] + 10.*t1077*t1078*t1079*var2[5] + 10.*t1073*t1074*t1075*var2[7] + 5.*t1067*t1068*t1071*var2[9] + t1063*t1065*var2[11];
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

  double *var1,*var2,*var3;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 3)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Three input(s) required (var1,var2,var3).");
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
    ( !(mrows == 12 && ncols == 1) && 
      !(mrows == 1 && ncols == 12))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
    ( !(mrows == 2 && ncols == 1) && 
      !(mrows == 1 && ncols == 2))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var3 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
  var3 = mxGetPr(prhs[2]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2,var3);


}

#else // MATLAB_MEX_FILE

#include "yd_VConstraint_Swing.hh"

namespace SymFunction
{

void yd_VConstraint_Swing_raw(double *p_output1, const double *var1,const double *var2,const double *var3)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3);

}

}

#endif // MATLAB_MEX_FILE
