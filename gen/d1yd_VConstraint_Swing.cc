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
static void output1(double *p_output1,const double *var1,const double *var2,const double *var3,const double *var4)
{
  double t1066;
  double t1072;
  double t1076;
  double t1080;
  double t1082;
  double t1086;
  double t1087;
  double t1088;
  double t1089;
  double t1090;
  double t1093;
  double t1094;
  double t1095;
  double t1098;
  double t1099;
  double t1102;
  double t1084;
  double t1085;
  double t1091;
  double t1092;
  double t1096;
  double t1097;
  double t1100;
  double t1101;
  double t1103;
  double t1104;
  double t1105;
  double t1109;
  double t1110;
  double t1111;
  double t1112;
  double t1113;
  double t1114;
  double t1115;
  double t1116;
  double t1117;
  double t1118;
  double t1119;
  t1066 = -1.*var4[0];
  t1072 = t1066 + var4[1];
  t1076 = Power(t1072,-5);
  t1080 = t1066 + var1[2] + var1[3];
  t1082 = Power(t1080,4);
  t1086 = Power(t1072,-4);
  t1087 = Power(t1080,3);
  t1088 = 1/t1072;
  t1089 = -1.*t1088*t1080;
  t1090 = 1. + t1089;
  t1093 = Power(t1072,-3);
  t1094 = Power(t1080,2);
  t1095 = Power(t1090,2);
  t1098 = Power(t1072,-2);
  t1099 = Power(t1090,3);
  t1102 = Power(t1090,4);
  t1084 = -5.*var3[8]*t1076*t1082;
  t1085 = 5.*var3[10]*t1076*t1082;
  t1091 = -20.*var3[6]*t1086*t1087*t1090;
  t1092 = 20.*var3[8]*t1086*t1087*t1090;
  t1096 = -30.*var3[4]*t1093*t1094*t1095;
  t1097 = 30.*var3[6]*t1093*t1094*t1095;
  t1100 = -20.*var3[2]*t1098*t1080*t1099;
  t1101 = 20.*var3[4]*t1098*t1080*t1099;
  t1103 = -5.*var3[0]*t1088*t1102;
  t1104 = 5.*var3[2]*t1088*t1102;
  t1105 = t1084 + t1085 + t1091 + t1092 + t1096 + t1097 + t1100 + t1101 + t1103 + t1104;
  t1109 = -5.*var3[9]*t1076*t1082;
  t1110 = 5.*var3[11]*t1076*t1082;
  t1111 = -20.*var3[7]*t1086*t1087*t1090;
  t1112 = 20.*var3[9]*t1086*t1087*t1090;
  t1113 = -30.*var3[5]*t1093*t1094*t1095;
  t1114 = 30.*var3[7]*t1093*t1094*t1095;
  t1115 = -20.*var3[3]*t1098*t1080*t1099;
  t1116 = 20.*var3[5]*t1098*t1080*t1099;
  t1117 = -5.*var3[1]*t1088*t1102;
  t1118 = 5.*var3[3]*t1088*t1102;
  t1119 = t1109 + t1110 + t1111 + t1112 + t1113 + t1114 + t1115 + t1116 + t1117 + t1118;
  p_output1[0]=t1105*var2[2] + t1105*var2[3];
  p_output1[1]=t1119*var2[2] + t1119*var2[3];
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

  double *var1,*var2,*var3,*var4;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 4)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Four input(s) required (var1,var2,var3,var4).");
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
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
    ( !(mrows == 12 && ncols == 1) && 
      !(mrows == 1 && ncols == 12))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var3 is wrong.");
    }
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
    ( !(mrows == 2 && ncols == 1) && 
      !(mrows == 1 && ncols == 2))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var4 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
  var3 = mxGetPr(prhs[2]);
  var4 = mxGetPr(prhs[3]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2,var3,var4);


}

#else // MATLAB_MEX_FILE

#include "d1yd_VConstraint_Swing.hh"

namespace SymFunction
{

void d1yd_VConstraint_Swing_raw(double *p_output1, const double *var1,const double *var2,const double *var3,const double *var4)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3, var4);

}

}

#endif // MATLAB_MEX_FILE
