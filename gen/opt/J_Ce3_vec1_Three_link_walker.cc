/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:51 GMT-04:00
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
  double t3637;
  double t8285;
  double t3638;
  double t8286;
  double t8319;
  double t8320;
  double t8321;
  double t8312;
  double t8313;
  double t8318;
  double t8277;
  double t8287;
  double t8288;
  double t8323;
  double t8324;
  double t8325;
  double t8301;
  double t8304;
  double t8332;
  double t8333;
  double t8334;
  double t8329;
  double t8330;
  double t8331;
  double t8303;
  double t8305;
  double t8306;
  double t8336;
  double t8337;
  double t8338;
  double t8289;
  double t8359;
  double t8360;
  double t8361;
  double t8362;
  double t8363;
  double t8322;
  double t8326;
  double t8327;
  double t8328;
  double t8344;
  double t8345;
  double t8346;
  double t8347;
  double t8348;
  double t8349;
  double t8311;
  double t8365;
  double t8366;
  double t8367;
  double t8368;
  double t8369;
  double t8335;
  double t8339;
  double t8340;
  double t8341;
  double t8350;
  double t8351;
  double t8352;
  double t8353;
  double t8354;
  double t8355;
  double t8397;
  double t8398;
  double t8399;
  double t8400;
  double t8401;
  double t8402;
  double t8403;
  double t8404;
  double t8405;
  double t8407;
  double t8408;
  double t8409;
  double t8410;
  double t8375;
  double t8376;
  double t8377;
  double t8378;
  double t8379;
  double t8380;
  double t8381;
  double t8389;
  double t8391;
  double t8392;
  double t8420;
  double t8415;
  double t8382;
  double t8383;
  double t8384;
  double t8385;
  double t8386;
  double t8387;
  double t8388;
  double t8390;
  double t8393;
  double t8394;
  double t8428;
  double t8416;
  t3637 = Cos(var1[2]);
  t8285 = Sin(var1[2]);
  t3638 = Cos(var1[3]);
  t8286 = Sin(var1[3]);
  t8319 = t3637*t3638;
  t8320 = -1.*t8285*t8286;
  t8321 = t8319 + t8320;
  t8312 = -1.*t3638*t8285;
  t8313 = -1.*t3637*t8286;
  t8318 = t8312 + t8313;
  t8277 = -1.*t3637*t3638;
  t8287 = t8285*t8286;
  t8288 = t8277 + t8287;
  t8323 = t3638*t8285;
  t8324 = t3637*t8286;
  t8325 = t8323 + t8324;
  t8301 = Cos(var1[4]);
  t8304 = Sin(var1[4]);
  t8332 = t3637*t8301;
  t8333 = -1.*t8285*t8304;
  t8334 = t8332 + t8333;
  t8329 = -1.*t8301*t8285;
  t8330 = -1.*t3637*t8304;
  t8331 = t8329 + t8330;
  t8303 = -1.*t3637*t8301;
  t8305 = t8285*t8304;
  t8306 = t8303 + t8305;
  t8336 = t8301*t8285;
  t8337 = t3637*t8304;
  t8338 = t8336 + t8337;
  t8289 = -1.25*var2[3]*t8288;
  t8359 = Power(t3638,2);
  t8360 = -0.5*t8359;
  t8361 = Power(t8286,2);
  t8362 = -0.5*t8361;
  t8363 = t8360 + t8362;
  t8322 = -15.*t8318*t8321;
  t8326 = -5.*t8325*t8321;
  t8327 = -15.*t8318*t8288;
  t8328 = -5.*t8325*t8288;
  t8344 = Power(t8318,2);
  t8345 = -10.*t8344;
  t8346 = -10.*t8318*t8325;
  t8347 = Power(t8321,2);
  t8348 = -10.*t8347;
  t8349 = -10.*t8321*t8288;
  t8311 = -1.25*var2[4]*t8306;
  t8365 = Power(t8301,2);
  t8366 = -0.5*t8365;
  t8367 = Power(t8304,2);
  t8368 = -0.5*t8367;
  t8369 = t8366 + t8368;
  t8335 = -15.*t8331*t8334;
  t8339 = -5.*t8338*t8334;
  t8340 = -15.*t8331*t8306;
  t8341 = -5.*t8338*t8306;
  t8350 = Power(t8331,2);
  t8351 = -10.*t8350;
  t8352 = -10.*t8331*t8338;
  t8353 = Power(t8334,2);
  t8354 = -10.*t8353;
  t8355 = -10.*t8334*t8306;
  t8397 = -5.*t8344;
  t8398 = -5.*t8318*t8325;
  t8399 = -5.*t8347;
  t8400 = -5.*t8321*t8288;
  t8401 = -5.*t8350;
  t8402 = -5.*t8331*t8338;
  t8403 = -5.*t8353;
  t8404 = -5.*t8334*t8306;
  t8405 = t8397 + t8398 + t8399 + t8400 + t8401 + t8402 + t8403 + t8404;
  t8407 = 10.*t8285;
  t8408 = -5.*t8318*t8363;
  t8409 = -5.*t8331*t8369;
  t8410 = t8407 + t8408 + t8409;
  t8375 = 2.5*var2[2]*t8288*t8363;
  t8376 = t8322 + t8326 + t8327 + t8328;
  t8377 = -0.5*var2[1]*t8376;
  t8378 = t8345 + t8346 + t8348 + t8349;
  t8379 = -0.5*var2[0]*t8378;
  t8380 = t8289 + t8375 + t8377 + t8379;
  t8381 = var2[0]*t8380;
  t8389 = -1.25*var2[3]*t8318;
  t8391 = -10.*t8318*t8321;
  t8392 = -10.*t8325*t8321;
  t8420 = t8397 + t8398 + t8399 + t8400;
  t8415 = -1.25*var2[0]*t8318;
  t8382 = 2.5*var2[2]*t8306*t8369;
  t8383 = t8335 + t8339 + t8340 + t8341;
  t8384 = -0.5*var2[1]*t8383;
  t8385 = t8351 + t8352 + t8354 + t8355;
  t8386 = -0.5*var2[0]*t8385;
  t8387 = t8311 + t8382 + t8384 + t8386;
  t8388 = var2[0]*t8387;
  t8390 = -1.25*var2[4]*t8331;
  t8393 = -10.*t8331*t8334;
  t8394 = -10.*t8338*t8334;
  t8428 = t8401 + t8402 + t8403 + t8404;
  t8416 = -1.25*var2[0]*t8331;
  p_output1[0]=var2[0]*(t8289 + t8311 - 0.5*(t8345 + t8346 + t8348 + t8349 + t8351 + t8352 + t8354 + t8355)*var2[0] - 0.5*(t8322 + t8326 + t8327 + t8328 + t8335 + t8339 + t8340 + t8341)*var2[1] - 0.5*(10.*t3637 - 5.*t8288*t8363 - 5.*t8306*t8369)*var2[2]);
  p_output1[1]=t8381;
  p_output1[2]=t8388;
  p_output1[3]=t8389 + t8390 - 1.*(t8391 + t8392 + t8393 + t8394)*var2[0] - 0.5*t8405*var2[1] - 0.5*t8410*var2[2];
  p_output1[4]=-0.5*t8405*var2[0];
  p_output1[5]=-0.5*t8410*var2[0];
  p_output1[6]=t8415;
  p_output1[7]=t8416;
  p_output1[8]=t8381;
  p_output1[9]=t8381;
  p_output1[10]=t8389 - 1.*(t8391 + t8392)*var2[0] - 0.5*t8420*var2[1] + 2.5*t8318*t8363*var2[2];
  p_output1[11]=-0.5*t8420*var2[0];
  p_output1[12]=2.5*t8318*t8363*var2[0];
  p_output1[13]=t8415;
  p_output1[14]=t8388;
  p_output1[15]=t8388;
  p_output1[16]=t8390 - 1.*(t8393 + t8394)*var2[0] - 0.5*t8428*var2[1] + 2.5*t8331*t8369*var2[2];
  p_output1[17]=-0.5*t8428*var2[0];
  p_output1[18]=2.5*t8331*t8369*var2[0];
  p_output1[19]=t8416;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 20, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_Ce3_vec1_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce3_vec1_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
