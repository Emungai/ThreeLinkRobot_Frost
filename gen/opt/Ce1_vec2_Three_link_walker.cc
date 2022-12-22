/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:34 GMT-04:00
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
  double t7282;
  double t7283;
  double t7301;
  double t7302;
  double t7303;
  double t7304;
  double t7305;
  double t7330;
  double t7331;
  double t7332;
  double t7351;
  double t7352;
  double t7353;
  double t7354;
  double t7355;
  double t7362;
  double t7363;
  double t7364;
  double t7306;
  double t7316;
  double t7317;
  double t7327;
  double t7328;
  double t7329;
  double t7338;
  double t7344;
  double t7345;
  double t7346;
  double t7347;
  double t7348;
  double t7356;
  double t7357;
  double t7358;
  double t7359;
  double t7360;
  double t7361;
  double t7365;
  double t7366;
  double t7367;
  double t7368;
  double t7369;
  double t7370;
  double t7377;
  double t7378;
  double t7381;
  double t7382;
  double t7389;
  double t7390;
  double t7391;
  double t7392;
  double t7393;
  double t7395;
  double t7396;
  double t7397;
  double t7398;
  double t7399;
  t7282 = Cos(var1[3]);
  t7283 = Sin(var1[2]);
  t7301 = -1.*t7282*t7283;
  t7302 = Cos(var1[2]);
  t7303 = Sin(var1[3]);
  t7304 = -1.*t7302*t7303;
  t7305 = t7301 + t7304;
  t7330 = t7302*t7282;
  t7331 = -1.*t7283*t7303;
  t7332 = t7330 + t7331;
  t7351 = Cos(var1[4]);
  t7352 = -1.*t7351*t7283;
  t7353 = Sin(var1[4]);
  t7354 = -1.*t7302*t7353;
  t7355 = t7352 + t7354;
  t7362 = t7302*t7351;
  t7363 = -1.*t7283*t7353;
  t7364 = t7362 + t7363;
  t7306 = Power(t7305,2);
  t7316 = 5.*t7306;
  t7317 = t7282*t7283;
  t7327 = t7302*t7303;
  t7328 = t7317 + t7327;
  t7329 = 5.*t7305*t7328;
  t7338 = Power(t7332,2);
  t7344 = 5.*t7338;
  t7345 = -1.*t7302*t7282;
  t7346 = t7283*t7303;
  t7347 = t7345 + t7346;
  t7348 = 5.*t7332*t7347;
  t7356 = Power(t7355,2);
  t7357 = 5.*t7356;
  t7358 = t7351*t7283;
  t7359 = t7302*t7353;
  t7360 = t7358 + t7359;
  t7361 = 5.*t7355*t7360;
  t7365 = Power(t7364,2);
  t7366 = 5.*t7365;
  t7367 = -1.*t7302*t7351;
  t7368 = t7283*t7353;
  t7369 = t7367 + t7368;
  t7370 = 5.*t7364*t7369;
  t7377 = 10.*t7305*t7332;
  t7378 = 10.*t7305*t7347;
  t7381 = 10.*t7355*t7364;
  t7382 = 10.*t7355*t7369;
  t7389 = Power(t7282,2);
  t7390 = -0.5*t7389;
  t7391 = Power(t7303,2);
  t7392 = -0.5*t7391;
  t7393 = t7390 + t7392;
  t7395 = Power(t7351,2);
  t7396 = -0.5*t7395;
  t7397 = Power(t7353,2);
  t7398 = -0.5*t7397;
  t7399 = t7396 + t7398;
  p_output1[0]=var2[1]*(-0.5*(t7316 + t7329 + t7344 + t7348 + t7357 + t7361 + t7366 + t7370)*var2[2] - 0.5*(t7316 + t7329 + t7344 + t7348)*var2[3] - 0.5*(t7357 + t7361 + t7366 + t7370)*var2[4]);
  p_output1[1]=var2[1]*(-0.5*(t7377 + t7378 + t7381 + t7382)*var2[2] - 0.5*(t7377 + t7378)*var2[3] - 0.5*(t7381 + t7382)*var2[4]);
  p_output1[2]=var2[1]*(-0.5*(-10.*t7302 + 5.*t7347*t7393 + 5.*t7369*t7399)*var2[2] - 2.5*t7347*t7393*var2[3] - 2.5*t7369*t7399*var2[4]);
  p_output1[3]=var2[1]*(1.25*t7347*var2[2] + 1.25*t7347*var2[3]);
  p_output1[4]=var2[1]*(1.25*t7369*var2[2] + 1.25*t7369*var2[4]);
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

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce1_vec2_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
