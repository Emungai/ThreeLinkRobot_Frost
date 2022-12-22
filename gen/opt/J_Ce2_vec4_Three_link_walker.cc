/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:50 GMT-04:00
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
  double t2731;
  double t2744;
  double t3048;
  double t3249;
  double t3251;
  double t4742;
  double t4743;
  double t8151;
  double t8152;
  double t8153;
  double t8147;
  double t8148;
  double t8149;
  double t8155;
  double t8156;
  double t8158;
  double t4747;
  double t4748;
  double t4749;
  double t4750;
  double t8117;
  double t8130;
  double t8131;
  double t8154;
  double t8159;
  double t8160;
  double t8161;
  double t8162;
  double t8168;
  double t8179;
  double t8180;
  double t8181;
  double t8182;
  double t8187;
  double t8188;
  double t8193;
  double t8198;
  double t8199;
  double t8204;
  double t8205;
  double t8206;
  double t8207;
  double t8209;
  double t8210;
  double t8211;
  double t8212;
  double t8213;
  double t8224;
  double t8225;
  double t8226;
  double t8229;
  double t8230;
  double t8231;
  double t8232;
  double t8233;
  double t8234;
  double t8214;
  double t8235;
  double t8236;
  double t8243;
  double t8244;
  double t8249;
  double t8250;
  double t8215;
  double t8238;
  double t8258;
  double t8259;
  double t8260;
  double t8261;
  t2731 = Cos(var1[2]);
  t2744 = Cos(var1[3]);
  t3048 = -1.*t2731*t2744;
  t3249 = Sin(var1[2]);
  t3251 = Sin(var1[3]);
  t4742 = t3249*t3251;
  t4743 = t3048 + t4742;
  t8151 = t2731*t2744;
  t8152 = -1.*t3249*t3251;
  t8153 = t8151 + t8152;
  t8147 = -1.*t2744*t3249;
  t8148 = -1.*t2731*t3251;
  t8149 = t8147 + t8148;
  t8155 = t2744*t3249;
  t8156 = t2731*t3251;
  t8158 = t8155 + t8156;
  t4747 = 1.25*var2[3]*t4743;
  t4748 = Power(t2744,2);
  t4749 = -0.5*t4748;
  t4750 = Power(t3251,2);
  t8117 = -0.5*t4750;
  t8130 = t4749 + t8117;
  t8131 = -2.5*var2[2]*t4743*t8130;
  t8154 = 15.*t8149*t8153;
  t8159 = 5.*t8158*t8153;
  t8160 = 15.*t8149*t4743;
  t8161 = 5.*t8158*t4743;
  t8162 = t8154 + t8159 + t8160 + t8161;
  t8168 = -0.5*var2[1]*t8162;
  t8179 = Power(t8149,2);
  t8180 = 10.*t8179;
  t8181 = 10.*t8149*t8158;
  t8182 = Power(t8153,2);
  t8187 = 10.*t8182;
  t8188 = 10.*t8153*t4743;
  t8193 = t8180 + t8181 + t8187 + t8188;
  t8198 = -0.5*var2[0]*t8193;
  t8199 = t4747 + t8131 + t8168 + t8198;
  t8204 = var2[3]*t8199;
  t8205 = 10.*t8149*t8153;
  t8206 = 10.*t8158*t8153;
  t8207 = t8205 + t8206;
  t8209 = 5.*t8179;
  t8210 = 5.*t8149*t8158;
  t8211 = 5.*t8182;
  t8212 = 5.*t8153*t4743;
  t8213 = t8209 + t8210 + t8211 + t8212;
  t8224 = 1.25*var2[3]*t8158;
  t8225 = -2.5*var2[2]*t8158*t8130;
  t8226 = -0.5*var2[0]*t8162;
  t8229 = Power(t4743,2);
  t8230 = 10.*t8229;
  t8231 = t8180 + t8181 + t8188 + t8230;
  t8232 = -0.5*var2[1]*t8231;
  t8233 = t8224 + t8225 + t8226 + t8232;
  t8234 = var2[3]*t8233;
  t8214 = -0.5*var2[3]*t8213;
  t8235 = 10.*t8149*t4743;
  t8236 = t8205 + t8235;
  t8243 = -2.5*var2[1]*t8158*t8130;
  t8244 = -2.5*var2[0]*t4743*t8130;
  t8249 = t8243 + t8244;
  t8250 = var2[3]*t8249;
  t8215 = -2.5*var2[3]*t8149*t8130;
  t8238 = -2.5*var2[3]*t4743*t8130;
  t8258 = 1.25*var2[1]*t8158;
  t8259 = 1.25*var2[0]*t4743;
  t8260 = t8258 + t8259;
  t8261 = var2[3]*t8260;
  p_output1[0]=t8204;
  p_output1[1]=t8204;
  p_output1[2]=-0.5*t8207*var2[3];
  p_output1[3]=t8214;
  p_output1[4]=t8215;
  p_output1[5]=-0.5*t8207*var2[0] - 0.5*t8213*var2[1] - 2.5*t8130*t8149*var2[2] + 2.5*t8149*var2[3];
  p_output1[6]=t8234;
  p_output1[7]=t8234;
  p_output1[8]=t8214;
  p_output1[9]=-0.5*t8236*var2[3];
  p_output1[10]=t8238;
  p_output1[11]=t8131 - 0.5*t8213*var2[0] - 0.5*t8236*var2[1] + 2.5*t4743*var2[3];
  p_output1[12]=t8250;
  p_output1[13]=t8250;
  p_output1[14]=t8215;
  p_output1[15]=t8238;
  p_output1[16]=-2.5*t8130*t8149*var2[0] - 2.5*t4743*t8130*var2[1];
  p_output1[17]=t8261;
  p_output1[18]=t8261;
  p_output1[19]=1.25*t8149*var2[3];
  p_output1[20]=t4747;
  p_output1[21]=1.25*t8149*var2[0] + 1.25*t4743*var2[1];
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

#include "J_Ce2_vec4_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce2_vec4_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
