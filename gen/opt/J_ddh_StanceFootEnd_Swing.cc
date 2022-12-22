/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:33:09 GMT-04:00
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
  double t6861;
  double t6845;
  double t6856;
  double t6862;
  double t6866;
  double t6873;
  double t6874;
  double t6876;
  double t6877;
  double t6878;
  double t6879;
  double t6880;
  double t6884;
  double t6860;
  double t6863;
  double t6864;
  double t6865;
  double t6875;
  double t6885;
  double t6886;
  double t6887;
  double t6894;
  double t6895;
  double t6896;
  double t6897;
  double t6898;
  double t6899;
  double t6905;
  double t6906;
  double t6907;
  double t6908;
  double t6909;
  double t6888;
  double t6889;
  double t6890;
  double t6913;
  double t6914;
  double t6915;
  t6861 = Cos(var1[2]);
  t6845 = Cos(var1[3]);
  t6856 = Sin(var1[2]);
  t6862 = Sin(var1[3]);
  t6866 = -1.*t6861*t6845;
  t6873 = t6856*t6862;
  t6874 = t6866 + t6873;
  t6876 = t6845*t6856;
  t6877 = t6861*t6862;
  t6878 = t6876 + t6877;
  t6879 = var2[1]*t6878;
  t6880 = var2[0]*t6874;
  t6884 = t6879 + t6880;
  t6860 = -1.*t6845*t6856;
  t6863 = -1.*t6861*t6862;
  t6864 = t6860 + t6863;
  t6865 = var3[0]*t6864;
  t6875 = var3[1]*t6874;
  t6885 = var2[2]*t6884;
  t6886 = var2[3]*t6884;
  t6887 = t6865 + t6875 + t6885 + t6886;
  t6894 = var2[0]*t6864;
  t6895 = var2[1]*t6874;
  t6896 = t6894 + t6895;
  t6897 = t6861*t6845;
  t6898 = -1.*t6856*t6862;
  t6899 = t6897 + t6898;
  t6905 = var3[1]*t6864;
  t6906 = var3[0]*t6899;
  t6907 = var2[2]*t6896;
  t6908 = var2[3]*t6896;
  t6909 = t6905 + t6906 + t6907 + t6908;
  t6888 = var2[2]*t6864;
  t6889 = var2[3]*t6864;
  t6890 = t6888 + t6889;
  t6913 = var2[1]*t6864;
  t6914 = var2[0]*t6899;
  t6915 = t6913 + t6914;
  p_output1[0]=t6887;
  p_output1[1]=t6887;
  p_output1[2]=t6890;
  p_output1[3]=t6874*var2[2] + t6874*var2[3];
  p_output1[4]=t6896;
  p_output1[5]=t6896;
  p_output1[6]=t6899;
  p_output1[7]=t6864;
  p_output1[8]=-1.*Power(t6845,2) - 1.*Power(t6862,2);
  p_output1[9]=-1.;
  p_output1[10]=t6909;
  p_output1[11]=t6909;
  p_output1[12]=t6899*var2[2] + t6899*var2[3];
  p_output1[13]=t6890;
  p_output1[14]=t6915;
  p_output1[15]=t6915;
  p_output1[16]=t6878;
  p_output1[17]=t6899;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 18, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2,var3);


}

#else // MATLAB_MEX_FILE

#include "J_ddh_StanceFootEnd_Swing.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_ddh_StanceFootEnd_Swing_raw(double *p_output1, const double *var1,const double *var2,const double *var3)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3);

}

}

#endif // MATLAB_MEX_FILE
