/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:30 GMT-04:00
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
static void output1(double *p_output1,const double *var1,const double *var2,const double *var3,const double *var4,const double *var5,const double *var6,const double *var7,const double *var8)
{
  double t7098;
  double t7099;
  double t7109;
  double t7122;
  double t7137;
  double t7139;
  double t7143;
  double t7145;
  double t7147;
  double t7144;
  double t7150;
  double t7152;
  double t7155;
  double t7161;
  double t7166;
  double t7167;
  double t7170;
  double t7171;
  double t7176;
  double t7179;
  double t7174;
  double t7175;
  double t7188;
  double t7189;
  double t7192;
  double t7193;
  double t7196;
  double t7197;
  t7098 = 4.*var5[0];
  t7099 = t7098 + var7[0] + var3[0];
  t7109 = -1. + var8[0];
  t7122 = 1/t7109;
  t7137 = -1.*var1[0];
  t7139 = t7137 + var1[1];
  t7143 = -0.333333333333333*t7122*t7139;
  t7145 = 4.*var5[1];
  t7147 = t7145 + var7[1] + var3[1];
  t7144 = -1.33333333333333*t7122*t7139;
  t7150 = 4.*var5[2];
  t7152 = t7150 + var7[2] + var3[2];
  t7155 = 4.*var5[3];
  t7161 = t7155 + var7[3] + var3[3];
  t7166 = 4.*var5[4];
  t7167 = t7166 + var7[4] + var3[4];
  t7170 = -1.*var7[0];
  t7171 = t7170 + var3[0];
  t7176 = -1.*var7[1];
  t7179 = t7176 + var3[1];
  t7174 = -0.25*t7122*t7139;
  t7175 = 0.25*t7122*t7139;
  t7188 = -1.*var7[2];
  t7189 = t7188 + var3[2];
  t7192 = -1.*var7[3];
  t7193 = t7192 + var3[3];
  t7196 = -1.*var7[4];
  t7197 = t7196 + var3[4];
  p_output1[0]=0.333333333333333*t7099*t7122;
  p_output1[1]=-0.333333333333333*t7099*t7122;
  p_output1[2]=-1.;
  p_output1[3]=t7143;
  p_output1[4]=t7144;
  p_output1[5]=1.;
  p_output1[6]=t7143;
  p_output1[7]=0.333333333333333*t7122*t7147;
  p_output1[8]=-0.333333333333333*t7122*t7147;
  p_output1[9]=-1.;
  p_output1[10]=t7143;
  p_output1[11]=t7144;
  p_output1[12]=1.;
  p_output1[13]=t7143;
  p_output1[14]=0.333333333333333*t7122*t7152;
  p_output1[15]=-0.333333333333333*t7122*t7152;
  p_output1[16]=-1.;
  p_output1[17]=t7143;
  p_output1[18]=t7144;
  p_output1[19]=1.;
  p_output1[20]=t7143;
  p_output1[21]=0.333333333333333*t7122*t7161;
  p_output1[22]=-0.333333333333333*t7122*t7161;
  p_output1[23]=-1.;
  p_output1[24]=t7143;
  p_output1[25]=t7144;
  p_output1[26]=1.;
  p_output1[27]=t7143;
  p_output1[28]=0.333333333333333*t7122*t7167;
  p_output1[29]=-0.333333333333333*t7122*t7167;
  p_output1[30]=-1.;
  p_output1[31]=t7143;
  p_output1[32]=t7144;
  p_output1[33]=1.;
  p_output1[34]=t7143;
  p_output1[35]=0.25*t7122*t7171;
  p_output1[36]=-0.25*t7122*t7171;
  p_output1[37]=-0.5;
  p_output1[38]=t7174;
  p_output1[39]=1.;
  p_output1[40]=-0.5;
  p_output1[41]=t7175;
  p_output1[42]=0.25*t7122*t7179;
  p_output1[43]=-0.25*t7122*t7179;
  p_output1[44]=-0.5;
  p_output1[45]=t7174;
  p_output1[46]=1.;
  p_output1[47]=-0.5;
  p_output1[48]=t7175;
  p_output1[49]=0.25*t7122*t7189;
  p_output1[50]=-0.25*t7122*t7189;
  p_output1[51]=-0.5;
  p_output1[52]=t7174;
  p_output1[53]=1.;
  p_output1[54]=-0.5;
  p_output1[55]=t7175;
  p_output1[56]=0.25*t7122*t7193;
  p_output1[57]=-0.25*t7122*t7193;
  p_output1[58]=-0.5;
  p_output1[59]=t7174;
  p_output1[60]=1.;
  p_output1[61]=-0.5;
  p_output1[62]=t7175;
  p_output1[63]=0.25*t7122*t7197;
  p_output1[64]=-0.25*t7122*t7197;
  p_output1[65]=-0.5;
  p_output1[66]=t7174;
  p_output1[67]=1.;
  p_output1[68]=-0.5;
  p_output1[69]=t7175;
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

  double *var1,*var2,*var3,*var4,*var5,*var6,*var7,*var8;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 8)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Eight input(s) required (var1,var2,var3,var4,var5,var6,var7,var8).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 2 && ncols == 1) && 
      !(mrows == 1 && ncols == 2))) 
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
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var3 is wrong.");
    }
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var4 is wrong.");
    }
  mrows = mxGetM(prhs[4]);
  ncols = mxGetN(prhs[4]);
  if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var5 is wrong.");
    }
  mrows = mxGetM(prhs[5]);
  ncols = mxGetN(prhs[5]);
  if( !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var6 is wrong.");
    }
  mrows = mxGetM(prhs[6]);
  ncols = mxGetN(prhs[6]);
  if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var7 is wrong.");
    }
  mrows = mxGetM(prhs[7]);
  ncols = mxGetN(prhs[7]);
  if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) ||
    ( !(mrows == 1 && ncols == 1) && 
      !(mrows == 1 && ncols == 1))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var8 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
  var3 = mxGetPr(prhs[2]);
  var4 = mxGetPr(prhs[3]);
  var5 = mxGetPr(prhs[4]);
  var6 = mxGetPr(prhs[5]);
  var7 = mxGetPr(prhs[6]);
  var8 = mxGetPr(prhs[7]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 70, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2,var3,var4,var5,var6,var7,var8);


}

#else // MATLAB_MEX_FILE

#include "J_hs_int_x.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_hs_int_x_raw(double *p_output1, const double *var1,const double *var2,const double *var3,const double *var4,const double *var5,const double *var6,const double *var7,const double *var8)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3, var4, var5, var6, var7, var8);

}

}

#endif // MATLAB_MEX_FILE
