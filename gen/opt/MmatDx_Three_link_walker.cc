/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:33 GMT-04:00
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
  double t4621;
  double t7200;
  double t7172;
  double t7201;
  double t7173;
  double t7205;
  double t7206;
  double t7210;
  double t7214;
  double t7213;
  double t7217;
  double t7218;
  double t7232;
  double t7235;
  double t7236;
  double t7242;
  double t7243;
  double t7244;
  double t7224;
  double t7227;
  double t7228;
  double t7238;
  double t7239;
  double t7240;
  double t7231;
  double t7237;
  double t7241;
  double t7245;
  double t7246;
  double t7248;
  double t7249;
  double t7250;
  double t7251;
  double t7254;
  double t7255;
  double t7258;
  double t7259;
  double t7263;
  double t7264;
  double t7265;
  double t7266;
  double t7267;
  double t7269;
  double t7270;
  double t7271;
  double t7272;
  double t7273;
  double t7287;
  double t7288;
  double t7289;
  double t7290;
  double t7262;
  double t7268;
  double t7274;
  double t7275;
  t4621 = Cos(var1[2]);
  t7200 = Sin(var1[2]);
  t7172 = Cos(var1[3]);
  t7201 = Sin(var1[3]);
  t7173 = t4621*t7172;
  t7205 = -1.*t7200*t7201;
  t7206 = t7173 + t7205;
  t7210 = Cos(var1[4]);
  t7214 = Sin(var1[4]);
  t7213 = t4621*t7210;
  t7217 = -1.*t7200*t7214;
  t7218 = t7213 + t7217;
  t7232 = t7172*t7200;
  t7235 = t4621*t7201;
  t7236 = t7232 + t7235;
  t7242 = t7210*t7200;
  t7243 = t4621*t7214;
  t7244 = t7242 + t7243;
  t7224 = -1.*t7172*t7200;
  t7227 = -1.*t4621*t7201;
  t7228 = t7224 + t7227;
  t7238 = -1.*t7210*t7200;
  t7239 = -1.*t4621*t7214;
  t7240 = t7238 + t7239;
  t7231 = -5.*t7228*t7206;
  t7237 = -5.*t7236*t7206;
  t7241 = -5.*t7240*t7218;
  t7245 = -5.*t7244*t7218;
  t7246 = t7231 + t7237 + t7241 + t7245;
  t7248 = Power(t4621,2);
  t7249 = -25.*t7248;
  t7250 = Power(t7200,2);
  t7251 = -25.*t7250;
  t7254 = Power(t7206,2);
  t7255 = -5.*t7254;
  t7258 = Power(t7218,2);
  t7259 = -5.*t7258;
  t7263 = Power(t7172,2);
  t7264 = -0.5*t7263;
  t7265 = Power(t7201,2);
  t7266 = -0.5*t7265;
  t7267 = t7264 + t7266;
  t7269 = Power(t7210,2);
  t7270 = -0.5*t7269;
  t7271 = Power(t7214,2);
  t7272 = -0.5*t7271;
  t7273 = t7270 + t7272;
  t7287 = 10.*t7200;
  t7288 = -5.*t7228*t7267;
  t7289 = -5.*t7240*t7273;
  t7290 = t7287 + t7288 + t7289;
  t7262 = -10.*t4621;
  t7268 = -5.*t7206*t7267;
  t7274 = -5.*t7218*t7273;
  t7275 = t7262 + t7268 + t7274;
  p_output1[0]=(-5.*Power(t7236,2) - 5.*Power(t7244,2) + t7249 + t7251 + t7255 + t7259)*var2[0] + t7246*var2[1] + t7275*var2[2] + 2.5*t7206*var2[3] + 2.5*t7218*var2[4];
  p_output1[1]=t7246*var2[0] + (-5.*Power(t7228,2) - 5.*Power(t7240,2) + t7249 + t7251 + t7255 + t7259)*var2[1] + t7290*var2[2] + 2.5*t7228*var2[3] + 2.5*t7240*var2[4];
  p_output1[2]=t7275*var2[0] + t7290*var2[1] + (-4.47 - 5.*Power(t7267,2) - 5.*Power(t7273,2))*var2[2] + 2.5*t7267*var2[3] + 2.5*t7273*var2[4];
  p_output1[3]=2.5*t7206*var2[0] + 2.5*t7228*var2[1] + 2.5*t7267*var2[2] - 1.25*var2[3];
  p_output1[4]=2.5*t7218*var2[0] + 2.5*t7240*var2[1] + 2.5*t7273*var2[2] - 1.25*var2[4];
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

#include "MmatDx_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void MmatDx_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
