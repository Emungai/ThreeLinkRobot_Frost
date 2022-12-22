/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:36 GMT-05:00
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
  double t239;
  double t245;
  double t253;
  double t261;
  double t262;
  double t263;
  double t283;
  double t305;
  double t306;
  double t307;
  double t313;
  double t314;
  double t315;
  double t324;
  double t325;
  double t326;
  double t299;
  double t300;
  double t301;
  double t302;
  double t303;
  double t312;
  double t319;
  double t320;
  double t321;
  double t322;
  double t323;
  double t327;
  double t328;
  t239 = Cos(var1[3]);
  t245 = Sin(var1[2]);
  t253 = -1.*t239*t245;
  t261 = Cos(var1[2]);
  t262 = Sin(var1[3]);
  t263 = -1.*t261*t262;
  t283 = t253 + t263;
  t305 = t261*t239;
  t306 = -1.*t245*t262;
  t307 = t305 + t306;
  t313 = t239*t245;
  t314 = t261*t262;
  t315 = t313 + t314;
  t324 = -1.*t261*t239;
  t325 = t245*t262;
  t326 = t324 + t325;
  t299 = Power(t239,2);
  t300 = -0.5*t299;
  t301 = Power(t262,2);
  t302 = -0.5*t301;
  t303 = t300 + t302;
  t312 = 10.*t283*t307;
  t319 = Power(t283,2);
  t320 = 5.*t319;
  t321 = 5.*t283*t315;
  t322 = Power(t307,2);
  t323 = 5.*t322;
  t327 = 5.*t307*t326;
  t328 = t320 + t321 + t323 + t327;
  p_output1[0]=var2[3]*(-0.5*(t312 + 10.*t307*t315)*var2[0] - 0.5*t328*var2[1] - 2.5*t283*t303*var2[2] + 1.25*t283*var2[3]);
  p_output1[1]=var2[3]*(-0.5*t328*var2[0] - 0.5*(t312 + 10.*t283*t326)*var2[1] - 2.5*t303*t326*var2[2] + 1.25*t326*var2[3]);
  p_output1[2]=(-2.5*t283*t303*var2[0] - 2.5*t303*t326*var2[1])*var2[3];
  p_output1[3]=(1.25*t283*var2[0] + 1.25*t326*var2[1])*var2[3];
  p_output1[4]=0;
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

#include "Ce2_vec4_Three_link_walker.hh"

namespace SymFunction
{

void Ce2_vec4_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
