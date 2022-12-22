/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:33 GMT-05:00
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
  double t230;
  double t232;
  double t229;
  double t231;
  double t236;
  double t237;
  double t238;
  double t246;
  double t247;
  double t248;
  double t240;
  double t241;
  double t242;
  double t243;
  double t244;
  double t254;
  double t255;
  double t256;
  double t250;
  double t251;
  double t252;
  double t258;
  double t259;
  double t260;
  double t269;
  double t270;
  double t271;
  double t278;
  double t279;
  double t280;
  double t249;
  double t257;
  double t264;
  double t265;
  double t266;
  double t267;
  double t268;
  double t272;
  double t273;
  double t274;
  double t275;
  double t276;
  double t277;
  double t281;
  double t282;
  double t285;
  double t286;
  double t287;
  double t288;
  double t289;
  double t291;
  double t292;
  double t293;
  double t294;
  double t295;
  double t284;
  double t290;
  double t296;
  double t297;
  double t308;
  double t309;
  double t310;
  double t311;
  t230 = Sin(var1[2]);
  t232 = Cos(var1[2]);
  t229 = Cos(var1[3]);
  t231 = -1.*t229*t230;
  t236 = Sin(var1[3]);
  t237 = -1.*t232*t236;
  t238 = t231 + t237;
  t246 = t232*t229;
  t247 = -1.*t230*t236;
  t248 = t246 + t247;
  t240 = Cos(var1[4]);
  t241 = -1.*t240*t230;
  t242 = Sin(var1[4]);
  t243 = -1.*t232*t242;
  t244 = t241 + t243;
  t254 = t232*t240;
  t255 = -1.*t230*t242;
  t256 = t254 + t255;
  t250 = t229*t230;
  t251 = t232*t236;
  t252 = t250 + t251;
  t258 = t240*t230;
  t259 = t232*t242;
  t260 = t258 + t259;
  t269 = -1.*t232*t229;
  t270 = t230*t236;
  t271 = t269 + t270;
  t278 = -1.*t232*t240;
  t279 = t230*t242;
  t280 = t278 + t279;
  t249 = 10.*t238*t248;
  t257 = 10.*t244*t256;
  t264 = Power(t238,2);
  t265 = 5.*t264;
  t266 = 5.*t238*t252;
  t267 = Power(t248,2);
  t268 = 5.*t267;
  t272 = 5.*t248*t271;
  t273 = Power(t244,2);
  t274 = 5.*t273;
  t275 = 5.*t244*t260;
  t276 = Power(t256,2);
  t277 = 5.*t276;
  t281 = 5.*t256*t280;
  t282 = t265 + t266 + t268 + t272 + t274 + t275 + t277 + t281;
  t285 = Power(t229,2);
  t286 = -0.5*t285;
  t287 = Power(t236,2);
  t288 = -0.5*t287;
  t289 = t286 + t288;
  t291 = Power(t240,2);
  t292 = -0.5*t291;
  t293 = Power(t242,2);
  t294 = -0.5*t293;
  t295 = t292 + t294;
  t284 = -10.*t230;
  t290 = 5.*t238*t289;
  t296 = 5.*t244*t295;
  t297 = t284 + t290 + t296;
  t308 = -10.*t232;
  t309 = 5.*t271*t289;
  t310 = 5.*t280*t295;
  t311 = t308 + t309 + t310;
  p_output1[0]=var2[2]*(-0.5*(t249 + 10.*t248*t252 + t257 + 10.*t256*t260)*var2[0] - 0.5*t282*var2[1] - 0.5*t297*var2[2] + 1.25*t238*var2[3] + 1.25*t244*var2[4]);
  p_output1[1]=var2[2]*(-0.5*t282*var2[0] - 0.5*(t249 + t257 + 10.*t238*t271 + 10.*t244*t280)*var2[1] - 0.5*t311*var2[2] + 1.25*t271*var2[3] + 1.25*t280*var2[4]);
  p_output1[2]=(-0.5*t297*var2[0] - 0.5*t311*var2[1])*var2[2];
  p_output1[3]=(1.25*t238*var2[0] + 1.25*t271*var2[1])*var2[2];
  p_output1[4]=(1.25*t244*var2[0] + 1.25*t280*var2[1])*var2[2];
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

namespace SymFunction
{

void Ce2_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
