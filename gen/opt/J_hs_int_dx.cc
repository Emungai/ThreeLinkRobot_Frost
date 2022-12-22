/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:32 GMT-04:00
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
  double t7124;
  double t7148;
  double t7165;
  double t7168;
  double t7186;
  double t7187;
  double t7191;
  double t7195;
  double t7198;
  double t7194;
  double t7202;
  double t7203;
  double t7207;
  double t7208;
  double t7211;
  double t7212;
  double t7215;
  double t7216;
  double t7221;
  double t7222;
  double t7219;
  double t7220;
  double t7225;
  double t7226;
  double t7229;
  double t7230;
  double t7233;
  double t7234;
  t7124 = 4.*var5[0];
  t7148 = var3[0] + t7124 + var7[0];
  t7165 = -1. + var8[0];
  t7168 = 1/t7165;
  t7186 = -1.*var1[0];
  t7187 = t7186 + var1[1];
  t7191 = -0.333333333333333*t7168*t7187;
  t7195 = 4.*var5[1];
  t7198 = var3[1] + t7195 + var7[1];
  t7194 = -1.33333333333333*t7168*t7187;
  t7202 = 4.*var5[2];
  t7203 = var3[2] + t7202 + var7[2];
  t7207 = 4.*var5[3];
  t7208 = var3[3] + t7207 + var7[3];
  t7211 = 4.*var5[4];
  t7212 = var3[4] + t7211 + var7[4];
  t7215 = -1.*var7[0];
  t7216 = var3[0] + t7215;
  t7221 = -1.*var7[1];
  t7222 = var3[1] + t7221;
  t7219 = -0.25*t7168*t7187;
  t7220 = 0.25*t7168*t7187;
  t7225 = -1.*var7[2];
  t7226 = var3[2] + t7225;
  t7229 = -1.*var7[3];
  t7230 = var3[3] + t7229;
  t7233 = -1.*var7[4];
  t7234 = var3[4] + t7233;
  p_output1[0]=0.333333333333333*t7148*t7168;
  p_output1[1]=-0.333333333333333*t7148*t7168;
  p_output1[2]=-1.;
  p_output1[3]=t7191;
  p_output1[4]=t7194;
  p_output1[5]=1.;
  p_output1[6]=t7191;
  p_output1[7]=0.333333333333333*t7168*t7198;
  p_output1[8]=-0.333333333333333*t7168*t7198;
  p_output1[9]=-1.;
  p_output1[10]=t7191;
  p_output1[11]=t7194;
  p_output1[12]=1.;
  p_output1[13]=t7191;
  p_output1[14]=0.333333333333333*t7168*t7203;
  p_output1[15]=-0.333333333333333*t7168*t7203;
  p_output1[16]=-1.;
  p_output1[17]=t7191;
  p_output1[18]=t7194;
  p_output1[19]=1.;
  p_output1[20]=t7191;
  p_output1[21]=0.333333333333333*t7168*t7208;
  p_output1[22]=-0.333333333333333*t7168*t7208;
  p_output1[23]=-1.;
  p_output1[24]=t7191;
  p_output1[25]=t7194;
  p_output1[26]=1.;
  p_output1[27]=t7191;
  p_output1[28]=0.333333333333333*t7168*t7212;
  p_output1[29]=-0.333333333333333*t7168*t7212;
  p_output1[30]=-1.;
  p_output1[31]=t7191;
  p_output1[32]=t7194;
  p_output1[33]=1.;
  p_output1[34]=t7191;
  p_output1[35]=0.25*t7168*t7216;
  p_output1[36]=-0.25*t7168*t7216;
  p_output1[37]=-0.5;
  p_output1[38]=t7219;
  p_output1[39]=1.;
  p_output1[40]=-0.5;
  p_output1[41]=t7220;
  p_output1[42]=0.25*t7168*t7222;
  p_output1[43]=-0.25*t7168*t7222;
  p_output1[44]=-0.5;
  p_output1[45]=t7219;
  p_output1[46]=1.;
  p_output1[47]=-0.5;
  p_output1[48]=t7220;
  p_output1[49]=0.25*t7168*t7226;
  p_output1[50]=-0.25*t7168*t7226;
  p_output1[51]=-0.5;
  p_output1[52]=t7219;
  p_output1[53]=1.;
  p_output1[54]=-0.5;
  p_output1[55]=t7220;
  p_output1[56]=0.25*t7168*t7230;
  p_output1[57]=-0.25*t7168*t7230;
  p_output1[58]=-0.5;
  p_output1[59]=t7219;
  p_output1[60]=1.;
  p_output1[61]=-0.5;
  p_output1[62]=t7220;
  p_output1[63]=0.25*t7168*t7234;
  p_output1[64]=-0.25*t7168*t7234;
  p_output1[65]=-0.5;
  p_output1[66]=t7219;
  p_output1[67]=1.;
  p_output1[68]=-0.5;
  p_output1[69]=t7220;
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

#include "J_hs_int_dx.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_hs_int_dx_raw(double *p_output1, const double *var1,const double *var2,const double *var3,const double *var4,const double *var5,const double *var6,const double *var7,const double *var8)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3, var4, var5, var6, var7, var8);

}

}

#endif // MATLAB_MEX_FILE
