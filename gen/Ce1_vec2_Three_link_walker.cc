/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:21 GMT-05:00
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
  double t87;
  double t88;
  double t102;
  double t103;
  double t104;
  double t105;
  double t106;
  double t131;
  double t132;
  double t133;
  double t152;
  double t153;
  double t154;
  double t155;
  double t156;
  double t163;
  double t164;
  double t165;
  double t107;
  double t117;
  double t118;
  double t128;
  double t129;
  double t130;
  double t139;
  double t145;
  double t146;
  double t147;
  double t148;
  double t149;
  double t157;
  double t158;
  double t159;
  double t160;
  double t161;
  double t162;
  double t166;
  double t167;
  double t168;
  double t169;
  double t170;
  double t171;
  double t178;
  double t179;
  double t182;
  double t183;
  double t190;
  double t191;
  double t192;
  double t193;
  double t194;
  double t196;
  double t197;
  double t198;
  double t199;
  double t200;
  t87 = Cos(var1[3]);
  t88 = Sin(var1[2]);
  t102 = -1.*t87*t88;
  t103 = Cos(var1[2]);
  t104 = Sin(var1[3]);
  t105 = -1.*t103*t104;
  t106 = t102 + t105;
  t131 = t103*t87;
  t132 = -1.*t88*t104;
  t133 = t131 + t132;
  t152 = Cos(var1[4]);
  t153 = -1.*t152*t88;
  t154 = Sin(var1[4]);
  t155 = -1.*t103*t154;
  t156 = t153 + t155;
  t163 = t103*t152;
  t164 = -1.*t88*t154;
  t165 = t163 + t164;
  t107 = Power(t106,2);
  t117 = 5.*t107;
  t118 = t87*t88;
  t128 = t103*t104;
  t129 = t118 + t128;
  t130 = 5.*t106*t129;
  t139 = Power(t133,2);
  t145 = 5.*t139;
  t146 = -1.*t103*t87;
  t147 = t88*t104;
  t148 = t146 + t147;
  t149 = 5.*t133*t148;
  t157 = Power(t156,2);
  t158 = 5.*t157;
  t159 = t152*t88;
  t160 = t103*t154;
  t161 = t159 + t160;
  t162 = 5.*t156*t161;
  t166 = Power(t165,2);
  t167 = 5.*t166;
  t168 = -1.*t103*t152;
  t169 = t88*t154;
  t170 = t168 + t169;
  t171 = 5.*t165*t170;
  t178 = 10.*t106*t133;
  t179 = 10.*t106*t148;
  t182 = 10.*t156*t165;
  t183 = 10.*t156*t170;
  t190 = Power(t87,2);
  t191 = -0.5*t190;
  t192 = Power(t104,2);
  t193 = -0.5*t192;
  t194 = t191 + t193;
  t196 = Power(t152,2);
  t197 = -0.5*t196;
  t198 = Power(t154,2);
  t199 = -0.5*t198;
  t200 = t197 + t199;
  p_output1[0]=var2[1]*(-0.5*(t117 + t130 + t145 + t149 + t158 + t162 + t167 + t171)*var2[2] - 0.5*(t117 + t130 + t145 + t149)*var2[3] - 0.5*(t158 + t162 + t167 + t171)*var2[4]);
  p_output1[1]=var2[1]*(-0.5*(t178 + t179 + t182 + t183)*var2[2] - 0.5*(t178 + t179)*var2[3] - 0.5*(t182 + t183)*var2[4]);
  p_output1[2]=var2[1]*(-0.5*(-10.*t103 + 5.*t148*t194 + 5.*t170*t200)*var2[2] - 2.5*t148*t194*var2[3] - 2.5*t170*t200*var2[4]);
  p_output1[3]=var2[1]*(1.25*t148*var2[2] + 1.25*t148*var2[3]);
  p_output1[4]=var2[1]*(1.25*t170*var2[2] + 1.25*t170*var2[4]);
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

#include "Ce1_vec2_Three_link_walker.hh"

namespace SymFunction
{

void Ce1_vec2_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
