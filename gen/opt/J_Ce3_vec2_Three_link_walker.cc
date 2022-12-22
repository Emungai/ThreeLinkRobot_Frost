/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:52 GMT-04:00
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
  double t3315;
  double t8343;
  double t3293;
  double t8356;
  double t8342;
  double t8357;
  double t8358;
  double t8412;
  double t8413;
  double t8414;
  double t8396;
  double t8406;
  double t8411;
  double t8419;
  double t8421;
  double t8422;
  double t8370;
  double t8372;
  double t8371;
  double t8373;
  double t8374;
  double t8429;
  double t8430;
  double t8431;
  double t8425;
  double t8426;
  double t8427;
  double t8434;
  double t8435;
  double t8436;
  double t8364;
  double t8456;
  double t8457;
  double t8458;
  double t8459;
  double t8460;
  double t8417;
  double t8418;
  double t8423;
  double t8424;
  double t8441;
  double t8442;
  double t8443;
  double t8444;
  double t8445;
  double t8446;
  double t8395;
  double t8462;
  double t8463;
  double t8464;
  double t8465;
  double t8466;
  double t8432;
  double t8433;
  double t8437;
  double t8438;
  double t8447;
  double t8448;
  double t8449;
  double t8450;
  double t8451;
  double t8452;
  double t8486;
  double t8487;
  double t8488;
  double t8489;
  double t8490;
  double t8491;
  double t8492;
  double t8493;
  double t8494;
  double t8495;
  double t8496;
  double t8507;
  double t8508;
  double t8509;
  double t8510;
  double t8472;
  double t8473;
  double t8474;
  double t8475;
  double t8476;
  double t8477;
  double t8478;
  double t8498;
  double t8500;
  double t8501;
  double t8516;
  double t8514;
  double t8479;
  double t8480;
  double t8481;
  double t8482;
  double t8483;
  double t8484;
  double t8485;
  double t8499;
  double t8502;
  double t8503;
  double t8524;
  double t8515;
  t3315 = Sin(var1[2]);
  t8343 = Cos(var1[2]);
  t3293 = Cos(var1[3]);
  t8356 = Sin(var1[3]);
  t8342 = t3293*t3315;
  t8357 = t8343*t8356;
  t8358 = t8342 + t8357;
  t8412 = t8343*t3293;
  t8413 = -1.*t3315*t8356;
  t8414 = t8412 + t8413;
  t8396 = -1.*t3293*t3315;
  t8406 = -1.*t8343*t8356;
  t8411 = t8396 + t8406;
  t8419 = -1.*t8343*t3293;
  t8421 = t3315*t8356;
  t8422 = t8419 + t8421;
  t8370 = Cos(var1[4]);
  t8372 = Sin(var1[4]);
  t8371 = t8370*t3315;
  t8373 = t8343*t8372;
  t8374 = t8371 + t8373;
  t8429 = t8343*t8370;
  t8430 = -1.*t3315*t8372;
  t8431 = t8429 + t8430;
  t8425 = -1.*t8370*t3315;
  t8426 = -1.*t8343*t8372;
  t8427 = t8425 + t8426;
  t8434 = -1.*t8343*t8370;
  t8435 = t3315*t8372;
  t8436 = t8434 + t8435;
  t8364 = -1.25*var2[3]*t8358;
  t8456 = Power(t3293,2);
  t8457 = -0.5*t8456;
  t8458 = Power(t8356,2);
  t8459 = -0.5*t8458;
  t8460 = t8457 + t8459;
  t8417 = -15.*t8411*t8414;
  t8418 = -5.*t8358*t8414;
  t8423 = -15.*t8411*t8422;
  t8424 = -5.*t8358*t8422;
  t8441 = Power(t8411,2);
  t8442 = -10.*t8441;
  t8443 = -10.*t8411*t8358;
  t8444 = -10.*t8414*t8422;
  t8445 = Power(t8422,2);
  t8446 = -10.*t8445;
  t8395 = -1.25*var2[4]*t8374;
  t8462 = Power(t8370,2);
  t8463 = -0.5*t8462;
  t8464 = Power(t8372,2);
  t8465 = -0.5*t8464;
  t8466 = t8463 + t8465;
  t8432 = -15.*t8427*t8431;
  t8433 = -5.*t8374*t8431;
  t8437 = -15.*t8427*t8436;
  t8438 = -5.*t8374*t8436;
  t8447 = Power(t8427,2);
  t8448 = -10.*t8447;
  t8449 = -10.*t8427*t8374;
  t8450 = -10.*t8431*t8436;
  t8451 = Power(t8436,2);
  t8452 = -10.*t8451;
  t8486 = -5.*t8441;
  t8487 = -5.*t8411*t8358;
  t8488 = Power(t8414,2);
  t8489 = -5.*t8488;
  t8490 = -5.*t8414*t8422;
  t8491 = -5.*t8447;
  t8492 = -5.*t8427*t8374;
  t8493 = Power(t8431,2);
  t8494 = -5.*t8493;
  t8495 = -5.*t8431*t8436;
  t8496 = t8486 + t8487 + t8489 + t8490 + t8491 + t8492 + t8494 + t8495;
  t8507 = 10.*t8343;
  t8508 = -5.*t8422*t8460;
  t8509 = -5.*t8436*t8466;
  t8510 = t8507 + t8508 + t8509;
  t8472 = 2.5*var2[2]*t8358*t8460;
  t8473 = t8417 + t8418 + t8423 + t8424;
  t8474 = -0.5*var2[0]*t8473;
  t8475 = t8442 + t8443 + t8444 + t8446;
  t8476 = -0.5*var2[1]*t8475;
  t8477 = t8364 + t8472 + t8474 + t8476;
  t8478 = var2[1]*t8477;
  t8498 = -1.25*var2[3]*t8422;
  t8500 = -10.*t8411*t8414;
  t8501 = -10.*t8411*t8422;
  t8516 = t8486 + t8487 + t8489 + t8490;
  t8514 = -1.25*var2[1]*t8422;
  t8479 = 2.5*var2[2]*t8374*t8466;
  t8480 = t8432 + t8433 + t8437 + t8438;
  t8481 = -0.5*var2[0]*t8480;
  t8482 = t8448 + t8449 + t8450 + t8452;
  t8483 = -0.5*var2[1]*t8482;
  t8484 = t8395 + t8479 + t8481 + t8483;
  t8485 = var2[1]*t8484;
  t8499 = -1.25*var2[4]*t8436;
  t8502 = -10.*t8427*t8431;
  t8503 = -10.*t8427*t8436;
  t8524 = t8491 + t8492 + t8494 + t8495;
  t8515 = -1.25*var2[1]*t8436;
  p_output1[0]=var2[1]*(t8364 + t8395 - 0.5*(t8417 + t8418 + t8423 + t8424 + t8432 + t8433 + t8437 + t8438)*var2[0] - 0.5*(t8442 + t8443 + t8444 + t8446 + t8448 + t8449 + t8450 + t8452)*var2[1] - 0.5*(-10.*t3315 - 5.*t8358*t8460 - 5.*t8374*t8466)*var2[2]);
  p_output1[1]=t8478;
  p_output1[2]=t8485;
  p_output1[3]=-0.5*t8496*var2[1];
  p_output1[4]=t8498 + t8499 - 0.5*t8496*var2[0] - 1.*(t8500 + t8501 + t8502 + t8503)*var2[1] - 0.5*t8510*var2[2];
  p_output1[5]=-0.5*t8510*var2[1];
  p_output1[6]=t8514;
  p_output1[7]=t8515;
  p_output1[8]=t8478;
  p_output1[9]=t8478;
  p_output1[10]=-0.5*t8516*var2[1];
  p_output1[11]=t8498 - 0.5*t8516*var2[0] - 1.*(t8500 + t8501)*var2[1] + 2.5*t8422*t8460*var2[2];
  p_output1[12]=2.5*t8422*t8460*var2[1];
  p_output1[13]=t8514;
  p_output1[14]=t8485;
  p_output1[15]=t8485;
  p_output1[16]=-0.5*t8524*var2[1];
  p_output1[17]=t8499 - 0.5*t8524*var2[0] - 1.*(t8502 + t8503)*var2[1] + 2.5*t8436*t8466*var2[2];
  p_output1[18]=2.5*t8436*t8466*var2[1];
  p_output1[19]=t8515;
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

#include "J_Ce3_vec2_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce3_vec2_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
