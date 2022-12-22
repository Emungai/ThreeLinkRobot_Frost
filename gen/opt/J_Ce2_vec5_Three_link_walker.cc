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
  double t3094;
  double t3099;
  double t4633;
  double t4642;
  double t4659;
  double t4667;
  double t4669;
  double t8222;
  double t8223;
  double t8237;
  double t8216;
  double t8217;
  double t8218;
  double t8240;
  double t8241;
  double t8242;
  double t4670;
  double t4671;
  double t4673;
  double t4674;
  double t4681;
  double t4682;
  double t8208;
  double t8239;
  double t8251;
  double t8256;
  double t8257;
  double t8262;
  double t8263;
  double t8264;
  double t8265;
  double t8266;
  double t8267;
  double t8268;
  double t8269;
  double t8270;
  double t8271;
  double t8272;
  double t8273;
  double t8274;
  double t8275;
  double t8276;
  double t8278;
  double t8279;
  double t8280;
  double t8281;
  double t8282;
  double t8290;
  double t8291;
  double t8292;
  double t8293;
  double t8294;
  double t8295;
  double t8296;
  double t8297;
  double t8298;
  double t8283;
  double t8299;
  double t8300;
  double t8307;
  double t8308;
  double t8309;
  double t8310;
  double t8284;
  double t8302;
  double t8314;
  double t8315;
  double t8316;
  double t8317;
  t3094 = Cos(var1[2]);
  t3099 = Cos(var1[4]);
  t4633 = -1.*t3094*t3099;
  t4642 = Sin(var1[2]);
  t4659 = Sin(var1[4]);
  t4667 = t4642*t4659;
  t4669 = t4633 + t4667;
  t8222 = t3094*t3099;
  t8223 = -1.*t4642*t4659;
  t8237 = t8222 + t8223;
  t8216 = -1.*t3099*t4642;
  t8217 = -1.*t3094*t4659;
  t8218 = t8216 + t8217;
  t8240 = t3099*t4642;
  t8241 = t3094*t4659;
  t8242 = t8240 + t8241;
  t4670 = 1.25*var2[4]*t4669;
  t4671 = Power(t3099,2);
  t4673 = -0.5*t4671;
  t4674 = Power(t4659,2);
  t4681 = -0.5*t4674;
  t4682 = t4673 + t4681;
  t8208 = -2.5*var2[2]*t4669*t4682;
  t8239 = 15.*t8218*t8237;
  t8251 = 5.*t8242*t8237;
  t8256 = 15.*t8218*t4669;
  t8257 = 5.*t8242*t4669;
  t8262 = t8239 + t8251 + t8256 + t8257;
  t8263 = -0.5*var2[1]*t8262;
  t8264 = Power(t8218,2);
  t8265 = 10.*t8264;
  t8266 = 10.*t8218*t8242;
  t8267 = Power(t8237,2);
  t8268 = 10.*t8267;
  t8269 = 10.*t8237*t4669;
  t8270 = t8265 + t8266 + t8268 + t8269;
  t8271 = -0.5*var2[0]*t8270;
  t8272 = t4670 + t8208 + t8263 + t8271;
  t8273 = var2[4]*t8272;
  t8274 = 10.*t8218*t8237;
  t8275 = 10.*t8242*t8237;
  t8276 = t8274 + t8275;
  t8278 = 5.*t8264;
  t8279 = 5.*t8218*t8242;
  t8280 = 5.*t8267;
  t8281 = 5.*t8237*t4669;
  t8282 = t8278 + t8279 + t8280 + t8281;
  t8290 = 1.25*var2[4]*t8242;
  t8291 = -2.5*var2[2]*t8242*t4682;
  t8292 = -0.5*var2[0]*t8262;
  t8293 = Power(t4669,2);
  t8294 = 10.*t8293;
  t8295 = t8265 + t8266 + t8269 + t8294;
  t8296 = -0.5*var2[1]*t8295;
  t8297 = t8290 + t8291 + t8292 + t8296;
  t8298 = var2[4]*t8297;
  t8283 = -0.5*var2[4]*t8282;
  t8299 = 10.*t8218*t4669;
  t8300 = t8274 + t8299;
  t8307 = -2.5*var2[1]*t8242*t4682;
  t8308 = -2.5*var2[0]*t4669*t4682;
  t8309 = t8307 + t8308;
  t8310 = var2[4]*t8309;
  t8284 = -2.5*var2[4]*t8218*t4682;
  t8302 = -2.5*var2[4]*t4669*t4682;
  t8314 = 1.25*var2[1]*t8242;
  t8315 = 1.25*var2[0]*t4669;
  t8316 = t8314 + t8315;
  t8317 = var2[4]*t8316;
  p_output1[0]=t8273;
  p_output1[1]=t8273;
  p_output1[2]=-0.5*t8276*var2[4];
  p_output1[3]=t8283;
  p_output1[4]=t8284;
  p_output1[5]=-0.5*t8276*var2[0] - 0.5*t8282*var2[1] - 2.5*t4682*t8218*var2[2] + 2.5*t8218*var2[4];
  p_output1[6]=t8298;
  p_output1[7]=t8298;
  p_output1[8]=t8283;
  p_output1[9]=-0.5*t8300*var2[4];
  p_output1[10]=t8302;
  p_output1[11]=t8208 - 0.5*t8282*var2[0] - 0.5*t8300*var2[1] + 2.5*t4669*var2[4];
  p_output1[12]=t8310;
  p_output1[13]=t8310;
  p_output1[14]=t8284;
  p_output1[15]=t8302;
  p_output1[16]=-2.5*t4682*t8218*var2[0] - 2.5*t4669*t4682*var2[1];
  p_output1[17]=t8317;
  p_output1[18]=t8317;
  p_output1[19]=1.25*t8218*var2[4];
  p_output1[20]=t4670;
  p_output1[21]=1.25*t8218*var2[0] + 1.25*t4669*var2[1];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 22, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_Ce2_vec5_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce2_vec5_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
