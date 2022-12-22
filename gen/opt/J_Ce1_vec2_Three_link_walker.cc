/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:43 GMT-04:00
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
  double t4174;
  double t4142;
  double t4161;
  double t4175;
  double t4198;
  double t4200;
  double t4205;
  double t4172;
  double t4180;
  double t4186;
  double t4207;
  double t4210;
  double t4219;
  double t4237;
  double t4243;
  double t4256;
  double t4647;
  double t4650;
  double t7874;
  double t7875;
  double t7876;
  double t4648;
  double t4651;
  double t7873;
  double t7886;
  double t7887;
  double t7888;
  double t7890;
  double t7894;
  double t7898;
  double t4206;
  double t4231;
  double t4257;
  double t4264;
  double t7885;
  double t7889;
  double t7900;
  double t7901;
  double t4645;
  double t4646;
  double t7902;
  double t7903;
  double t7926;
  double t7932;
  double t7938;
  double t7940;
  double t7941;
  double t7942;
  double t7957;
  double t7958;
  double t7959;
  double t7960;
  double t7961;
  double t7962;
  double t7965;
  double t7943;
  double t7963;
  double t7971;
  double t7976;
  double t7977;
  double t7978;
  double t7979;
  double t7986;
  double t7987;
  double t7993;
  double t7994;
  double t7995;
  double t7980;
  double t7985;
  double t7997;
  double t7998;
  double t8009;
  double t8010;
  double t8013;
  double t8014;
  double t8017;
  double t8011;
  double t8015;
  double t8023;
  double t8024;
  double t8025;
  double t8026;
  double t8027;
  double t8029;
  double t8030;
  double t8031;
  double t8032;
  double t8033;
  double t8028;
  double t8034;
  double t8050;
  double t8051;
  double t8052;
  double t8053;
  double t8059;
  double t8060;
  double t8061;
  double t8062;
  double t8066;
  double t8067;
  double t8068;
  double t8069;
  double t8070;
  double t8074;
  t4174 = Cos(var1[2]);
  t4142 = Cos(var1[3]);
  t4161 = Sin(var1[2]);
  t4175 = Sin(var1[3]);
  t4198 = t4174*t4142;
  t4200 = -1.*t4161*t4175;
  t4205 = t4198 + t4200;
  t4172 = -1.*t4142*t4161;
  t4180 = -1.*t4174*t4175;
  t4186 = t4172 + t4180;
  t4207 = t4142*t4161;
  t4210 = t4174*t4175;
  t4219 = t4207 + t4210;
  t4237 = -1.*t4174*t4142;
  t4243 = t4161*t4175;
  t4256 = t4237 + t4243;
  t4647 = Cos(var1[4]);
  t4650 = Sin(var1[4]);
  t7874 = t4174*t4647;
  t7875 = -1.*t4161*t4650;
  t7876 = t7874 + t7875;
  t4648 = -1.*t4647*t4161;
  t4651 = -1.*t4174*t4650;
  t7873 = t4648 + t4651;
  t7886 = t4647*t4161;
  t7887 = t4174*t4650;
  t7888 = t7886 + t7887;
  t7890 = -1.*t4174*t4647;
  t7894 = t4161*t4650;
  t7898 = t7890 + t7894;
  t4206 = 15.*t4186*t4205;
  t4231 = 5.*t4219*t4205;
  t4257 = 15.*t4186*t4256;
  t4264 = 5.*t4219*t4256;
  t7885 = 15.*t7873*t7876;
  t7889 = 5.*t7888*t7876;
  t7900 = 15.*t7873*t7898;
  t7901 = 5.*t7888*t7898;
  t4645 = t4206 + t4231 + t4257 + t4264;
  t4646 = -0.5*var2[3]*t4645;
  t7902 = t7885 + t7889 + t7900 + t7901;
  t7903 = -0.5*var2[4]*t7902;
  t7926 = Power(t4186,2);
  t7932 = 5.*t7926;
  t7938 = 5.*t4186*t4219;
  t7940 = Power(t4205,2);
  t7941 = 5.*t7940;
  t7942 = 5.*t4205*t4256;
  t7957 = Power(t7873,2);
  t7958 = 5.*t7957;
  t7959 = 5.*t7873*t7888;
  t7960 = Power(t7876,2);
  t7961 = 5.*t7960;
  t7962 = 5.*t7876*t7898;
  t7965 = t7932 + t7938 + t7941 + t7942 + t7958 + t7959 + t7961 + t7962;
  t7943 = t7932 + t7938 + t7941 + t7942;
  t7963 = t7958 + t7959 + t7961 + t7962;
  t7971 = 10.*t7926;
  t7976 = 10.*t4186*t4219;
  t7977 = 10.*t4205*t4256;
  t7978 = Power(t4256,2);
  t7979 = 10.*t7978;
  t7986 = 10.*t7957;
  t7987 = 10.*t7873*t7888;
  t7993 = 10.*t7876*t7898;
  t7994 = Power(t7898,2);
  t7995 = 10.*t7994;
  t7980 = t7971 + t7976 + t7977 + t7979;
  t7985 = -0.5*var2[3]*t7980;
  t7997 = t7986 + t7987 + t7993 + t7995;
  t7998 = -0.5*var2[4]*t7997;
  t8009 = 10.*t4186*t4205;
  t8010 = 10.*t4186*t4256;
  t8013 = 10.*t7873*t7876;
  t8014 = 10.*t7873*t7898;
  t8017 = t8009 + t8010 + t8013 + t8014;
  t8011 = t8009 + t8010;
  t8015 = t8013 + t8014;
  t8023 = Power(t4142,2);
  t8024 = -0.5*t8023;
  t8025 = Power(t4175,2);
  t8026 = -0.5*t8025;
  t8027 = t8024 + t8026;
  t8029 = Power(t4647,2);
  t8030 = -0.5*t8029;
  t8031 = Power(t4650,2);
  t8032 = -0.5*t8031;
  t8033 = t8030 + t8032;
  t8028 = -2.5*var2[3]*t4219*t8027;
  t8034 = -2.5*var2[4]*t7888*t8033;
  t8050 = -10.*t4174;
  t8051 = 5.*t4256*t8027;
  t8052 = 5.*t7898*t8033;
  t8053 = t8050 + t8051 + t8052;
  t8059 = 1.25*var2[2]*t4219;
  t8060 = 1.25*var2[3]*t4219;
  t8061 = t8059 + t8060;
  t8062 = var2[1]*t8061;
  t8066 = 1.25*var2[1]*t4256;
  t8067 = 1.25*var2[2]*t7888;
  t8068 = 1.25*var2[4]*t7888;
  t8069 = t8067 + t8068;
  t8070 = var2[1]*t8069;
  t8074 = 1.25*var2[1]*t7898;
  p_output1[0]=var2[1]*(t4646 + t7903 - 0.5*(t4206 + t4231 + t4257 + t4264 + t7885 + t7889 + t7900 + t7901)*var2[2]);
  p_output1[1]=var2[1]*(t4646 - 0.5*t4645*var2[2]);
  p_output1[2]=var2[1]*(t7903 - 0.5*t7902*var2[2]);
  p_output1[3]=-0.5*t7965*var2[2] - 0.5*t7943*var2[3] - 0.5*t7963*var2[4];
  p_output1[4]=-0.5*t7965*var2[1];
  p_output1[5]=-0.5*t7943*var2[1];
  p_output1[6]=-0.5*t7963*var2[1];
  p_output1[7]=var2[1]*(t7985 + t7998 - 0.5*(t7971 + t7976 + t7977 + t7979 + t7986 + t7987 + t7993 + t7995)*var2[2]);
  p_output1[8]=var2[1]*(t7985 - 0.5*t7980*var2[2]);
  p_output1[9]=var2[1]*(t7998 - 0.5*t7997*var2[2]);
  p_output1[10]=-0.5*t8017*var2[2] - 0.5*t8011*var2[3] - 0.5*t8015*var2[4];
  p_output1[11]=-0.5*t8017*var2[1];
  p_output1[12]=-0.5*t8011*var2[1];
  p_output1[13]=-0.5*t8015*var2[1];
  p_output1[14]=var2[1]*(t8028 + t8034 - 0.5*(10.*t4161 + 5.*t4219*t8027 + 5.*t7888*t8033)*var2[2]);
  p_output1[15]=var2[1]*(t8028 - 2.5*t4219*t8027*var2[2]);
  p_output1[16]=var2[1]*(t8034 - 2.5*t7888*t8033*var2[2]);
  p_output1[17]=-0.5*t8053*var2[2] - 2.5*t4256*t8027*var2[3] - 2.5*t7898*t8033*var2[4];
  p_output1[18]=-0.5*t8053*var2[1];
  p_output1[19]=-2.5*t4256*t8027*var2[1];
  p_output1[20]=-2.5*t7898*t8033*var2[1];
  p_output1[21]=t8062;
  p_output1[22]=t8062;
  p_output1[23]=1.25*t4256*var2[2] + 1.25*t4256*var2[3];
  p_output1[24]=t8066;
  p_output1[25]=t8066;
  p_output1[26]=t8070;
  p_output1[27]=t8070;
  p_output1[28]=1.25*t7898*var2[2] + 1.25*t7898*var2[4];
  p_output1[29]=t8074;
  p_output1[30]=t8074;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 31, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_Ce1_vec2_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce1_vec2_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
