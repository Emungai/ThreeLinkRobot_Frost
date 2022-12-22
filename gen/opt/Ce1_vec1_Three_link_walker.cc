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
  double t7252;
  double t7209;
  double t7223;
  double t7253;
  double t7260;
  double t7261;
  double t7276;
  double t7284;
  double t7286;
  double t7293;
  double t7294;
  double t7295;
  double t7247;
  double t7256;
  double t7257;
  double t7277;
  double t7278;
  double t7279;
  double t7280;
  double t7281;
  double t7285;
  double t7291;
  double t7292;
  double t7296;
  double t7297;
  double t7298;
  double t7299;
  double t7300;
  double t7307;
  double t7308;
  double t7309;
  double t7310;
  double t7311;
  double t7312;
  double t7313;
  double t7314;
  double t7315;
  double t7318;
  double t7319;
  double t7320;
  double t7321;
  double t7322;
  double t7323;
  double t7324;
  double t7325;
  double t7326;
  double t7333;
  double t7334;
  double t7335;
  double t7336;
  double t7337;
  double t7339;
  double t7340;
  double t7341;
  double t7342;
  double t7343;
  t7252 = Cos(var1[2]);
  t7209 = Cos(var1[3]);
  t7223 = Sin(var1[2]);
  t7253 = Sin(var1[3]);
  t7260 = t7252*t7209;
  t7261 = -1.*t7223*t7253;
  t7276 = t7260 + t7261;
  t7284 = Cos(var1[4]);
  t7286 = Sin(var1[4]);
  t7293 = t7252*t7284;
  t7294 = -1.*t7223*t7286;
  t7295 = t7293 + t7294;
  t7247 = -1.*t7209*t7223;
  t7256 = -1.*t7252*t7253;
  t7257 = t7247 + t7256;
  t7277 = 10.*t7257*t7276;
  t7278 = t7209*t7223;
  t7279 = t7252*t7253;
  t7280 = t7278 + t7279;
  t7281 = 10.*t7280*t7276;
  t7285 = -1.*t7284*t7223;
  t7291 = -1.*t7252*t7286;
  t7292 = t7285 + t7291;
  t7296 = 10.*t7292*t7295;
  t7297 = t7284*t7223;
  t7298 = t7252*t7286;
  t7299 = t7297 + t7298;
  t7300 = 10.*t7299*t7295;
  t7307 = Power(t7257,2);
  t7308 = 5.*t7307;
  t7309 = 5.*t7257*t7280;
  t7310 = Power(t7276,2);
  t7311 = 5.*t7310;
  t7312 = -1.*t7252*t7209;
  t7313 = t7223*t7253;
  t7314 = t7312 + t7313;
  t7315 = 5.*t7276*t7314;
  t7318 = Power(t7292,2);
  t7319 = 5.*t7318;
  t7320 = 5.*t7292*t7299;
  t7321 = Power(t7295,2);
  t7322 = 5.*t7321;
  t7323 = -1.*t7252*t7284;
  t7324 = t7223*t7286;
  t7325 = t7323 + t7324;
  t7326 = 5.*t7295*t7325;
  t7333 = Power(t7209,2);
  t7334 = -0.5*t7333;
  t7335 = Power(t7253,2);
  t7336 = -0.5*t7335;
  t7337 = t7334 + t7336;
  t7339 = Power(t7284,2);
  t7340 = -0.5*t7339;
  t7341 = Power(t7286,2);
  t7342 = -0.5*t7341;
  t7343 = t7340 + t7342;
  p_output1[0]=var2[0]*(-0.5*(t7277 + t7281 + t7296 + t7300)*var2[2] - 0.5*(t7277 + t7281)*var2[3] - 0.5*(t7296 + t7300)*var2[4]);
  p_output1[1]=var2[0]*(-0.5*(t7308 + t7309 + t7311 + t7315 + t7319 + t7320 + t7322 + t7326)*var2[2] - 0.5*(t7308 + t7309 + t7311 + t7315)*var2[3] - 0.5*(t7319 + t7320 + t7322 + t7326)*var2[4]);
  p_output1[2]=var2[0]*(-0.5*(-10.*t7223 + 5.*t7257*t7337 + 5.*t7292*t7343)*var2[2] - 2.5*t7257*t7337*var2[3] - 2.5*t7292*t7343*var2[4]);
  p_output1[3]=var2[0]*(1.25*t7257*var2[2] + 1.25*t7257*var2[3]);
  p_output1[4]=var2[0]*(1.25*t7292*var2[2] + 1.25*t7292*var2[4]);
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

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce1_vec1_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
