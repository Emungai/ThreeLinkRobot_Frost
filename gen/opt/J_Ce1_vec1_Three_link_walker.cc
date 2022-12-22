/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:42 GMT-04:00
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
  double t7767;
  double t7768;
  double t7788;
  double t7803;
  double t7804;
  double t7805;
  double t7806;
  double t7815;
  double t7816;
  double t7821;
  double t7849;
  double t7850;
  double t7851;
  double t7852;
  double t7853;
  double t7860;
  double t7863;
  double t7864;
  double t7807;
  double t7809;
  double t7810;
  double t7811;
  double t7812;
  double t7813;
  double t7822;
  double t7823;
  double t7824;
  double t7825;
  double t7841;
  double t7842;
  double t7854;
  double t7855;
  double t7856;
  double t7857;
  double t7858;
  double t7859;
  double t7865;
  double t7866;
  double t7867;
  double t7868;
  double t7869;
  double t7870;
  double t7843;
  double t7848;
  double t7871;
  double t7872;
  double t7891;
  double t7892;
  double t7895;
  double t7896;
  double t7899;
  double t7893;
  double t7897;
  double t7905;
  double t7906;
  double t7907;
  double t7908;
  double t7911;
  double t7912;
  double t7913;
  double t7914;
  double t7909;
  double t7910;
  double t7915;
  double t7916;
  double t7927;
  double t7928;
  double t7929;
  double t7930;
  double t7933;
  double t7934;
  double t7935;
  double t7936;
  double t7939;
  double t7931;
  double t7937;
  double t7945;
  double t7946;
  double t7947;
  double t7948;
  double t7949;
  double t7951;
  double t7952;
  double t7953;
  double t7954;
  double t7955;
  double t7950;
  double t7956;
  double t7972;
  double t7973;
  double t7974;
  double t7975;
  double t7981;
  double t7982;
  double t7983;
  double t7984;
  double t7988;
  double t7989;
  double t7990;
  double t7991;
  double t7992;
  double t7996;
  t7767 = Cos(var1[3]);
  t7768 = Sin(var1[2]);
  t7788 = -1.*t7767*t7768;
  t7803 = Cos(var1[2]);
  t7804 = Sin(var1[3]);
  t7805 = -1.*t7803*t7804;
  t7806 = t7788 + t7805;
  t7815 = t7803*t7767;
  t7816 = -1.*t7768*t7804;
  t7821 = t7815 + t7816;
  t7849 = Cos(var1[4]);
  t7850 = -1.*t7849*t7768;
  t7851 = Sin(var1[4]);
  t7852 = -1.*t7803*t7851;
  t7853 = t7850 + t7852;
  t7860 = t7803*t7849;
  t7863 = -1.*t7768*t7851;
  t7864 = t7860 + t7863;
  t7807 = Power(t7806,2);
  t7809 = 10.*t7807;
  t7810 = t7767*t7768;
  t7811 = t7803*t7804;
  t7812 = t7810 + t7811;
  t7813 = 10.*t7806*t7812;
  t7822 = Power(t7821,2);
  t7823 = 10.*t7822;
  t7824 = -1.*t7803*t7767;
  t7825 = t7768*t7804;
  t7841 = t7824 + t7825;
  t7842 = 10.*t7821*t7841;
  t7854 = Power(t7853,2);
  t7855 = 10.*t7854;
  t7856 = t7849*t7768;
  t7857 = t7803*t7851;
  t7858 = t7856 + t7857;
  t7859 = 10.*t7853*t7858;
  t7865 = Power(t7864,2);
  t7866 = 10.*t7865;
  t7867 = -1.*t7803*t7849;
  t7868 = t7768*t7851;
  t7869 = t7867 + t7868;
  t7870 = 10.*t7864*t7869;
  t7843 = t7809 + t7813 + t7823 + t7842;
  t7848 = -0.5*var2[3]*t7843;
  t7871 = t7855 + t7859 + t7866 + t7870;
  t7872 = -0.5*var2[4]*t7871;
  t7891 = 10.*t7806*t7821;
  t7892 = 10.*t7812*t7821;
  t7895 = 10.*t7853*t7864;
  t7896 = 10.*t7858*t7864;
  t7899 = t7891 + t7892 + t7895 + t7896;
  t7893 = t7891 + t7892;
  t7897 = t7895 + t7896;
  t7905 = 15.*t7806*t7821;
  t7906 = 5.*t7812*t7821;
  t7907 = 15.*t7806*t7841;
  t7908 = 5.*t7812*t7841;
  t7911 = 15.*t7853*t7864;
  t7912 = 5.*t7858*t7864;
  t7913 = 15.*t7853*t7869;
  t7914 = 5.*t7858*t7869;
  t7909 = t7905 + t7906 + t7907 + t7908;
  t7910 = -0.5*var2[3]*t7909;
  t7915 = t7911 + t7912 + t7913 + t7914;
  t7916 = -0.5*var2[4]*t7915;
  t7927 = 5.*t7807;
  t7928 = 5.*t7806*t7812;
  t7929 = 5.*t7822;
  t7930 = 5.*t7821*t7841;
  t7933 = 5.*t7854;
  t7934 = 5.*t7853*t7858;
  t7935 = 5.*t7865;
  t7936 = 5.*t7864*t7869;
  t7939 = t7927 + t7928 + t7929 + t7930 + t7933 + t7934 + t7935 + t7936;
  t7931 = t7927 + t7928 + t7929 + t7930;
  t7937 = t7933 + t7934 + t7935 + t7936;
  t7945 = Power(t7767,2);
  t7946 = -0.5*t7945;
  t7947 = Power(t7804,2);
  t7948 = -0.5*t7947;
  t7949 = t7946 + t7948;
  t7951 = Power(t7849,2);
  t7952 = -0.5*t7951;
  t7953 = Power(t7851,2);
  t7954 = -0.5*t7953;
  t7955 = t7952 + t7954;
  t7950 = -2.5*var2[3]*t7841*t7949;
  t7956 = -2.5*var2[4]*t7869*t7955;
  t7972 = -10.*t7768;
  t7973 = 5.*t7806*t7949;
  t7974 = 5.*t7853*t7955;
  t7975 = t7972 + t7973 + t7974;
  t7981 = 1.25*var2[2]*t7841;
  t7982 = 1.25*var2[3]*t7841;
  t7983 = t7981 + t7982;
  t7984 = var2[0]*t7983;
  t7988 = 1.25*var2[0]*t7806;
  t7989 = 1.25*var2[2]*t7869;
  t7990 = 1.25*var2[4]*t7869;
  t7991 = t7989 + t7990;
  t7992 = var2[0]*t7991;
  t7996 = 1.25*var2[0]*t7853;
  p_output1[0]=var2[0]*(t7848 + t7872 - 0.5*(t7809 + t7813 + t7823 + t7842 + t7855 + t7859 + t7866 + t7870)*var2[2]);
  p_output1[1]=var2[0]*(t7848 - 0.5*t7843*var2[2]);
  p_output1[2]=var2[0]*(t7872 - 0.5*t7871*var2[2]);
  p_output1[3]=-0.5*t7899*var2[2] - 0.5*t7893*var2[3] - 0.5*t7897*var2[4];
  p_output1[4]=-0.5*t7899*var2[0];
  p_output1[5]=-0.5*t7893*var2[0];
  p_output1[6]=-0.5*t7897*var2[0];
  p_output1[7]=var2[0]*(t7910 + t7916 - 0.5*(t7905 + t7906 + t7907 + t7908 + t7911 + t7912 + t7913 + t7914)*var2[2]);
  p_output1[8]=var2[0]*(t7910 - 0.5*t7909*var2[2]);
  p_output1[9]=var2[0]*(t7916 - 0.5*t7915*var2[2]);
  p_output1[10]=-0.5*t7939*var2[2] - 0.5*t7931*var2[3] - 0.5*t7937*var2[4];
  p_output1[11]=-0.5*t7939*var2[0];
  p_output1[12]=-0.5*t7931*var2[0];
  p_output1[13]=-0.5*t7937*var2[0];
  p_output1[14]=var2[0]*(t7950 + t7956 - 0.5*(-10.*t7803 + 5.*t7841*t7949 + 5.*t7869*t7955)*var2[2]);
  p_output1[15]=var2[0]*(t7950 - 2.5*t7841*t7949*var2[2]);
  p_output1[16]=var2[0]*(t7956 - 2.5*t7869*t7955*var2[2]);
  p_output1[17]=-0.5*t7975*var2[2] - 2.5*t7806*t7949*var2[3] - 2.5*t7853*t7955*var2[4];
  p_output1[18]=-0.5*t7975*var2[0];
  p_output1[19]=-2.5*t7806*t7949*var2[0];
  p_output1[20]=-2.5*t7853*t7955*var2[0];
  p_output1[21]=t7984;
  p_output1[22]=t7984;
  p_output1[23]=1.25*t7806*var2[2] + 1.25*t7806*var2[3];
  p_output1[24]=t7988;
  p_output1[25]=t7988;
  p_output1[26]=t7992;
  p_output1[27]=t7992;
  p_output1[28]=1.25*t7853*var2[2] + 1.25*t7853*var2[4];
  p_output1[29]=t7996;
  p_output1[30]=t7996;
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

#include "J_Ce1_vec1_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce1_vec1_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
