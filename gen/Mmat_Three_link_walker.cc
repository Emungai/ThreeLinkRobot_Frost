/*
 * Automatically Generated from Mathematica.
 * Wed 15 Jul 2020 15:09:02 GMT-04:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
#include<math.h>
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

inline double Sec(double x) { return 1/cos(x); }
inline double Csc(double x) { return 1/sin(x); }

#endif

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t6;
  double t3;
  double t9;
  double t11;
  double t21;
  double t23;
  double t16;
  double t17;
  double t18;
  double t10;
  double t12;
  double t13;
  double t28;
  double t29;
  double t30;
  double t22;
  double t24;
  double t25;
  double t34;
  double t35;
  double t36;
  double t37;
  double t38;
  double t39;
  double t40;
  double t41;
  double t42;
  double t43;
  double t44;
  double t4;
  double t5;
  double t7;
  double t8;
  double t19;
  double t20;
  double t31;
  double t32;
  double t46;
  double t47;
  double t48;
  double t49;
  double t50;
  double t52;
  double t53;
  double t54;
  double t55;
  double t56;
  double t45;
  double t51;
  double t57;
  double t58;
  double t66;
  double t67;
  double t68;
  double t69;
  double t59;
  double t70;
  double t77;
  double t60;
  double t71;
  double t78;
  t6 = Sin(var1[2]);
  t3 = Cos(var1[2]);
  t9 = Cos(var1[3]);
  t11 = Sin(var1[3]);
  t21 = Cos(var1[4]);
  t23 = Sin(var1[4]);
  t16 = t3*t9;
  t17 = -1.*t6*t11;
  t18 = t16 + t17;
  t10 = t9*t6;
  t12 = t3*t11;
  t13 = t10 + t12;
  t28 = t3*t21;
  t29 = -1.*t6*t23;
  t30 = t28 + t29;
  t22 = t21*t6;
  t24 = t3*t23;
  t25 = t22 + t24;
  t34 = -1.*t9*t6;
  t35 = -1.*t3*t11;
  t36 = t34 + t35;
  t37 = 5.*t36*t18;
  t38 = 5.*t13*t18;
  t39 = -1.*t21*t6;
  t40 = -1.*t3*t23;
  t41 = t39 + t40;
  t42 = 5.*t41*t30;
  t43 = 5.*t25*t30;
  t44 = t37 + t38 + t42 + t43;
  t4 = Power(t3,2);
  t5 = 25.*t4;
  t7 = Power(t6,2);
  t8 = 25.*t7;
  t19 = Power(t18,2);
  t20 = 5.*t19;
  t31 = Power(t30,2);
  t32 = 5.*t31;
  t46 = Power(t9,2);
  t47 = -0.5*t46;
  t48 = Power(t11,2);
  t49 = -0.5*t48;
  t50 = t47 + t49;
  t52 = Power(t21,2);
  t53 = -0.5*t52;
  t54 = Power(t23,2);
  t55 = -0.5*t54;
  t56 = t53 + t55;
  t45 = 10.*t3;
  t51 = 5.*t18*t50;
  t57 = 5.*t30*t56;
  t58 = t45 + t51 + t57;
  t66 = -10.*t6;
  t67 = 5.*t36*t50;
  t68 = 5.*t41*t56;
  t69 = t66 + t67 + t68;
  t59 = -2.5*t18;
  t70 = -2.5*t36;
  t77 = -2.5*t50;
  t60 = -2.5*t30;
  t71 = -2.5*t41;
  t78 = -2.5*t56;
  p_output1[0]=5.*Power(t13,2) + t20 + 5.*Power(t25,2) + t32 + t5 + t8;
  p_output1[1]=t44;
  p_output1[2]=t58;
  p_output1[3]=t59;
  p_output1[4]=t60;
  p_output1[5]=t44;
  p_output1[6]=t20 + t32 + 5.*Power(t36,2) + 5.*Power(t41,2) + t5 + t8;
  p_output1[7]=t69;
  p_output1[8]=t70;
  p_output1[9]=t71;
  p_output1[10]=t58;
  p_output1[11]=t69;
  p_output1[12]=4.47 + 5.*Power(t50,2) + 5.*Power(t56,2);
  p_output1[13]=t77;
  p_output1[14]=t78;
  p_output1[15]=t59;
  p_output1[16]=t70;
  p_output1[17]=t77;
  p_output1[18]=1.25;
  p_output1[19]=0;
  p_output1[20]=t60;
  p_output1[21]=t71;
  p_output1[22]=t78;
  p_output1[23]=0;
  p_output1[24]=1.25;
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

  double *var1;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "One input(s) required (var1).");
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

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 5, (mwSize) 5, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "Mmat_Three_link_walker.hh"

namespace SymFunction
{

void Mmat_Three_link_walker_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
