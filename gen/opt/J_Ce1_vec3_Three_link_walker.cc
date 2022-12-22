/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:44 GMT-04:00
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
  double t3252;
  double t3276;
  double t2747;
  double t3275;
  double t3316;
  double t3338;
  double t3265;
  double t3278;
  double t3279;
  double t3294;
  double t3300;
  double t3302;
  double t3306;
  double t3307;
  double t3337;
  double t3339;
  double t3340;
  double t3346;
  double t3626;
  double t3627;
  double t4263;
  double t7904;
  double t3310;
  double t7917;
  double t7969;
  double t7970;
  double t7999;
  double t8001;
  double t8002;
  double t8003;
  double t8020;
  double t8021;
  double t8022;
  double t8036;
  double t8037;
  double t8038;
  double t8035;
  double t8039;
  double t7918;
  double t7919;
  double t7920;
  double t7921;
  double t7925;
  double t7966;
  t3252 = Cos(var1[3]);
  t3276 = Sin(var1[3]);
  t2747 = Cos(var1[2]);
  t3275 = Sin(var1[2]);
  t3316 = Cos(var1[4]);
  t3338 = Sin(var1[4]);
  t3265 = -1.*t2747*t3252;
  t3278 = t3275*t3276;
  t3279 = t3265 + t3278;
  t3294 = Power(t3252,2);
  t3300 = -0.5*t3294;
  t3302 = Power(t3276,2);
  t3306 = -0.5*t3302;
  t3307 = t3300 + t3306;
  t3337 = -1.*t2747*t3316;
  t3339 = t3275*t3338;
  t3340 = t3337 + t3339;
  t3346 = Power(t3316,2);
  t3626 = -0.5*t3346;
  t3627 = Power(t3338,2);
  t4263 = -0.5*t3627;
  t7904 = t3626 + t4263;
  t3310 = -2.5*var2[3]*t3279*t3307;
  t7917 = -2.5*var2[4]*t3340*t7904;
  t7969 = -1.*t3252*t3275;
  t7970 = -1.*t2747*t3276;
  t7999 = t7969 + t7970;
  t8001 = -1.*t3316*t3275;
  t8002 = -1.*t2747*t3338;
  t8003 = t8001 + t8002;
  t8020 = t3252*t3275;
  t8021 = t2747*t3276;
  t8022 = t8020 + t8021;
  t8036 = t3316*t3275;
  t8037 = t2747*t3338;
  t8038 = t8036 + t8037;
  t8035 = -2.5*var2[3]*t8022*t3307;
  t8039 = -2.5*var2[4]*t8038*t7904;
  t7918 = -10.*t2747;
  t7919 = 5.*t3279*t3307;
  t7920 = 5.*t3340*t7904;
  t7921 = t7918 + t7919 + t7920;
  t7925 = -2.5*var2[2]*t3279*t3307;
  t7966 = -2.5*var2[2]*t3340*t7904;
  p_output1[0]=var2[2]*(t3310 + t7917 - 0.5*t7921*var2[2]);
  p_output1[1]=(t3310 + t7925)*var2[2];
  p_output1[2]=(t7917 + t7966)*var2[2];
  p_output1[3]=-1.*(-10.*t3275 + 5.*t3307*t7999 + 5.*t7904*t8003)*var2[2] - 2.5*t3307*t7999*var2[3] - 2.5*t7904*t8003*var2[4];
  p_output1[4]=-2.5*t3307*t7999*var2[2];
  p_output1[5]=-2.5*t7904*t8003*var2[2];
  p_output1[6]=var2[2]*(t8035 + t8039 - 0.5*(10.*t3275 + 5.*t3307*t8022 + 5.*t7904*t8038)*var2[2]);
  p_output1[7]=var2[2]*(t8035 - 2.5*t3307*t8022*var2[2]);
  p_output1[8]=var2[2]*(t8039 - 2.5*t7904*t8038*var2[2]);
  p_output1[9]=t3310 + t7917 - 1.*t7921*var2[2];
  p_output1[10]=t7925;
  p_output1[11]=t7966;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 12, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_Ce1_vec3_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce1_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
