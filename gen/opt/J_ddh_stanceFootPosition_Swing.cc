/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:45:02 GMT-04:00
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
static void output1(double *p_output1,const double *var1,const double *var2,const double *var3)
{
  double t8661;
  double t8662;
  double t8669;
  double t8670;
  double t8671;
  double t8672;
  double t8673;
  double t8682;
  double t8683;
  double t8684;
  double t8685;
  double t8686;
  double t8687;
  double t8680;
  double t8681;
  double t8688;
  double t8689;
  double t8690;
  double t8691;
  double t8692;
  double t8693;
  double t8694;
  double t8695;
  double t8696;
  double t8699;
  double t8700;
  double t8701;
  double t8702;
  double t8703;
  double t8704;
  double t8697;
  double t8698;
  double t8705;
  double t8706;
  double t8707;
  double t8708;
  double t8709;
  double t8710;
  t8661 = Cos(var1[3]);
  t8662 = Sin(var1[2]);
  t8669 = t8661*t8662;
  t8670 = Cos(var1[2]);
  t8671 = Sin(var1[3]);
  t8672 = t8670*t8671;
  t8673 = t8669 + t8672;
  t8682 = t8670*t8661;
  t8683 = -1.*t8662*t8671;
  t8684 = t8682 + t8683;
  t8685 = var2[2]*t8684;
  t8686 = var2[3]*t8684;
  t8687 = t8685 + t8686;
  t8680 = var3[2]*t8673;
  t8681 = var3[3]*t8673;
  t8688 = var2[2]*t8687;
  t8689 = var2[3]*t8687;
  t8690 = t8680 + t8681 + t8688 + t8689;
  t8691 = 2.*var2[2]*t8673;
  t8692 = 2.*var2[3]*t8673;
  t8693 = t8691 + t8692;
  t8694 = -1.*t8670*t8661;
  t8695 = t8662*t8671;
  t8696 = t8694 + t8695;
  t8699 = -1.*t8661*t8662;
  t8700 = -1.*t8670*t8671;
  t8701 = t8699 + t8700;
  t8702 = var2[2]*t8701;
  t8703 = var2[3]*t8701;
  t8704 = t8702 + t8703;
  t8697 = var3[2]*t8684;
  t8698 = var3[3]*t8684;
  t8705 = var2[2]*t8704;
  t8706 = var2[3]*t8704;
  t8707 = t8697 + t8698 + t8705 + t8706;
  t8708 = 2.*var2[2]*t8684;
  t8709 = 2.*var2[3]*t8684;
  t8710 = t8708 + t8709;
  p_output1[0]=t8690;
  p_output1[1]=t8690;
  p_output1[2]=t8693;
  p_output1[3]=t8693;
  p_output1[4]=1.;
  p_output1[5]=t8696;
  p_output1[6]=t8696;
  p_output1[7]=t8707;
  p_output1[8]=t8707;
  p_output1[9]=t8710;
  p_output1[10]=t8710;
  p_output1[11]=1.;
  p_output1[12]=t8673;
  p_output1[13]=t8673;
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

  double *var1,*var2,*var3;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 3)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Three input(s) required (var1,var2,var3).");
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
  mrows = mxGetM(prhs[2]);
  ncols = mxGetN(prhs[2]);
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
    ( !(mrows == 5 && ncols == 1) && 
      !(mrows == 1 && ncols == 5))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var3 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
  var3 = mxGetPr(prhs[2]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 14, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2,var3);


}

#else // MATLAB_MEX_FILE

#include "J_ddh_stanceFootPosition_Swing.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_ddh_stanceFootPosition_Swing_raw(double *p_output1, const double *var1,const double *var2,const double *var3)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3);

}

}

#endif // MATLAB_MEX_FILE
