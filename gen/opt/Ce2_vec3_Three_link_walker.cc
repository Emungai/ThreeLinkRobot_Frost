/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:36 GMT-04:00
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
  double t7429;
  double t7431;
  double t7428;
  double t7430;
  double t7435;
  double t7436;
  double t7437;
  double t7445;
  double t7446;
  double t7447;
  double t7439;
  double t7440;
  double t7441;
  double t7442;
  double t7443;
  double t7453;
  double t7454;
  double t7455;
  double t7449;
  double t7450;
  double t7451;
  double t7457;
  double t7458;
  double t7459;
  double t7468;
  double t7469;
  double t7470;
  double t7477;
  double t7478;
  double t7479;
  double t7448;
  double t7456;
  double t7463;
  double t7464;
  double t7465;
  double t7466;
  double t7467;
  double t7471;
  double t7472;
  double t7473;
  double t7474;
  double t7475;
  double t7476;
  double t7480;
  double t7481;
  double t7484;
  double t7485;
  double t7486;
  double t7487;
  double t7488;
  double t7490;
  double t7491;
  double t7492;
  double t7493;
  double t7494;
  double t7483;
  double t7489;
  double t7495;
  double t7496;
  double t7507;
  double t7508;
  double t7509;
  double t7510;
  t7429 = Sin(var1[2]);
  t7431 = Cos(var1[2]);
  t7428 = Cos(var1[3]);
  t7430 = -1.*t7428*t7429;
  t7435 = Sin(var1[3]);
  t7436 = -1.*t7431*t7435;
  t7437 = t7430 + t7436;
  t7445 = t7431*t7428;
  t7446 = -1.*t7429*t7435;
  t7447 = t7445 + t7446;
  t7439 = Cos(var1[4]);
  t7440 = -1.*t7439*t7429;
  t7441 = Sin(var1[4]);
  t7442 = -1.*t7431*t7441;
  t7443 = t7440 + t7442;
  t7453 = t7431*t7439;
  t7454 = -1.*t7429*t7441;
  t7455 = t7453 + t7454;
  t7449 = t7428*t7429;
  t7450 = t7431*t7435;
  t7451 = t7449 + t7450;
  t7457 = t7439*t7429;
  t7458 = t7431*t7441;
  t7459 = t7457 + t7458;
  t7468 = -1.*t7431*t7428;
  t7469 = t7429*t7435;
  t7470 = t7468 + t7469;
  t7477 = -1.*t7431*t7439;
  t7478 = t7429*t7441;
  t7479 = t7477 + t7478;
  t7448 = 10.*t7437*t7447;
  t7456 = 10.*t7443*t7455;
  t7463 = Power(t7437,2);
  t7464 = 5.*t7463;
  t7465 = 5.*t7437*t7451;
  t7466 = Power(t7447,2);
  t7467 = 5.*t7466;
  t7471 = 5.*t7447*t7470;
  t7472 = Power(t7443,2);
  t7473 = 5.*t7472;
  t7474 = 5.*t7443*t7459;
  t7475 = Power(t7455,2);
  t7476 = 5.*t7475;
  t7480 = 5.*t7455*t7479;
  t7481 = t7464 + t7465 + t7467 + t7471 + t7473 + t7474 + t7476 + t7480;
  t7484 = Power(t7428,2);
  t7485 = -0.5*t7484;
  t7486 = Power(t7435,2);
  t7487 = -0.5*t7486;
  t7488 = t7485 + t7487;
  t7490 = Power(t7439,2);
  t7491 = -0.5*t7490;
  t7492 = Power(t7441,2);
  t7493 = -0.5*t7492;
  t7494 = t7491 + t7493;
  t7483 = -10.*t7429;
  t7489 = 5.*t7437*t7488;
  t7495 = 5.*t7443*t7494;
  t7496 = t7483 + t7489 + t7495;
  t7507 = -10.*t7431;
  t7508 = 5.*t7470*t7488;
  t7509 = 5.*t7479*t7494;
  t7510 = t7507 + t7508 + t7509;
  p_output1[0]=var2[2]*(-0.5*(t7448 + 10.*t7447*t7451 + t7456 + 10.*t7455*t7459)*var2[0] - 0.5*t7481*var2[1] - 0.5*t7496*var2[2] + 1.25*t7437*var2[3] + 1.25*t7443*var2[4]);
  p_output1[1]=var2[2]*(-0.5*t7481*var2[0] - 0.5*(t7448 + t7456 + 10.*t7437*t7470 + 10.*t7443*t7479)*var2[1] - 0.5*t7510*var2[2] + 1.25*t7470*var2[3] + 1.25*t7479*var2[4]);
  p_output1[2]=(-0.5*t7496*var2[0] - 0.5*t7510*var2[1])*var2[2];
  p_output1[3]=(1.25*t7437*var2[0] + 1.25*t7470*var2[1])*var2[2];
  p_output1[4]=(1.25*t7443*var2[0] + 1.25*t7479*var2[1])*var2[2];
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

#include "Ce2_vec3_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void Ce2_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
