/*
 * Automatically Generated from Mathematica.
 * Tue 10 Dec 2019 01:08:23 GMT-05:00
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
  double t150;
  double t174;
  double t151;
  double t173;
  double t187;
  double t189;
  double t172;
  double t175;
  double t176;
  double t177;
  double t180;
  double t181;
  double t184;
  double t185;
  double t188;
  double t195;
  double t201;
  double t202;
  double t203;
  double t204;
  double t205;
  double t206;
  double t215;
  double t216;
  double t217;
  double t219;
  double t220;
  double t221;
  t150 = Cos(var1[3]);
  t174 = Sin(var1[3]);
  t151 = Sin(var1[2]);
  t173 = Cos(var1[2]);
  t187 = Cos(var1[4]);
  t189 = Sin(var1[4]);
  t172 = -1.*t150*t151;
  t175 = -1.*t173*t174;
  t176 = t172 + t175;
  t177 = Power(t150,2);
  t180 = -0.5*t177;
  t181 = Power(t174,2);
  t184 = -0.5*t181;
  t185 = t180 + t184;
  t188 = -1.*t187*t151;
  t195 = -1.*t173*t189;
  t201 = t188 + t195;
  t202 = Power(t187,2);
  t203 = -0.5*t202;
  t204 = Power(t189,2);
  t205 = -0.5*t204;
  t206 = t203 + t205;
  t215 = -1.*t173*t150;
  t216 = t151*t174;
  t217 = t215 + t216;
  t219 = -1.*t173*t187;
  t220 = t151*t189;
  t221 = t219 + t220;
  p_output1[0]=var2[2]*(-0.5*(-10.*t151 + 5.*t176*t185 + 5.*t201*t206)*var2[2] - 2.5*t176*t185*var2[3] - 2.5*t201*t206*var2[4]);
  p_output1[1]=var2[2]*(-0.5*(-10.*t173 + 5.*t185*t217 + 5.*t206*t221)*var2[2] - 2.5*t185*t217*var2[3] - 2.5*t206*t221*var2[4]);
  p_output1[2]=0;
  p_output1[3]=0;
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

#include "Ce1_vec3_Three_link_walker.hh"

namespace SymFunction
{

void Ce1_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
