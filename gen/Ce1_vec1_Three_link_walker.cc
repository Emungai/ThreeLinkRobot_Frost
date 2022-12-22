/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:19 GMT-05:00
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
  double t35;
  double t22;
  double t23;
  double t41;
  double t71;
  double t72;
  double t73;
  double t89;
  double t91;
  double t94;
  double t95;
  double t96;
  double t34;
  double t69;
  double t70;
  double t80;
  double t81;
  double t82;
  double t83;
  double t84;
  double t90;
  double t92;
  double t93;
  double t97;
  double t98;
  double t99;
  double t100;
  double t101;
  double t108;
  double t109;
  double t110;
  double t111;
  double t112;
  double t113;
  double t114;
  double t115;
  double t116;
  double t119;
  double t120;
  double t121;
  double t122;
  double t123;
  double t124;
  double t125;
  double t126;
  double t127;
  double t134;
  double t135;
  double t136;
  double t137;
  double t138;
  double t140;
  double t141;
  double t142;
  double t143;
  double t144;
  t35 = Cos(var1[2]);
  t22 = Cos(var1[3]);
  t23 = Sin(var1[2]);
  t41 = Sin(var1[3]);
  t71 = t35*t22;
  t72 = -1.*t23*t41;
  t73 = t71 + t72;
  t89 = Cos(var1[4]);
  t91 = Sin(var1[4]);
  t94 = t35*t89;
  t95 = -1.*t23*t91;
  t96 = t94 + t95;
  t34 = -1.*t22*t23;
  t69 = -1.*t35*t41;
  t70 = t34 + t69;
  t80 = 10.*t70*t73;
  t81 = t22*t23;
  t82 = t35*t41;
  t83 = t81 + t82;
  t84 = 10.*t83*t73;
  t90 = -1.*t89*t23;
  t92 = -1.*t35*t91;
  t93 = t90 + t92;
  t97 = 10.*t93*t96;
  t98 = t89*t23;
  t99 = t35*t91;
  t100 = t98 + t99;
  t101 = 10.*t100*t96;
  t108 = Power(t70,2);
  t109 = 5.*t108;
  t110 = 5.*t70*t83;
  t111 = Power(t73,2);
  t112 = 5.*t111;
  t113 = -1.*t35*t22;
  t114 = t23*t41;
  t115 = t113 + t114;
  t116 = 5.*t73*t115;
  t119 = Power(t93,2);
  t120 = 5.*t119;
  t121 = 5.*t93*t100;
  t122 = Power(t96,2);
  t123 = 5.*t122;
  t124 = -1.*t35*t89;
  t125 = t23*t91;
  t126 = t124 + t125;
  t127 = 5.*t96*t126;
  t134 = Power(t22,2);
  t135 = -0.5*t134;
  t136 = Power(t41,2);
  t137 = -0.5*t136;
  t138 = t135 + t137;
  t140 = Power(t89,2);
  t141 = -0.5*t140;
  t142 = Power(t91,2);
  t143 = -0.5*t142;
  t144 = t141 + t143;
  p_output1[0]=var2[0]*(-0.5*(t101 + t80 + t84 + t97)*var2[2] - 0.5*(t80 + t84)*var2[3] - 0.5*(t101 + t97)*var2[4]);
  p_output1[1]=var2[0]*(-0.5*(t109 + t110 + t112 + t116 + t120 + t121 + t123 + t127)*var2[2] - 0.5*(t109 + t110 + t112 + t116)*var2[3] - 0.5*(t120 + t121 + t123 + t127)*var2[4]);
  p_output1[2]=var2[0]*(-0.5*(-10.*t23 + 5.*t138*t70 + 5.*t144*t93)*var2[2] - 2.5*t138*t70*var2[3] - 2.5*t144*t93*var2[4]);
  p_output1[3]=var2[0]*(1.25*t70*var2[2] + 1.25*t70*var2[3]);
  p_output1[4]=var2[0]*(1.25*t93*var2[2] + 1.25*t93*var2[4]);
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

#include "Ce1_vec1_Three_link_walker.hh"

namespace SymFunction
{

void Ce1_vec1_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
