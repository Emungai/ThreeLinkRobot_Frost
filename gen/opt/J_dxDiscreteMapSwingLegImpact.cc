/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:45:15 GMT-04:00
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
static void output1(double *p_output1,const double *var1,const double *var2,const double *var3,const double *var4)
{
  double t2831;
  double t2820;
  double t2840;
  double t2842;
  double t2853;
  double t2855;
  double t2848;
  double t2849;
  double t2850;
  double t2841;
  double t2844;
  double t2845;
  double t2861;
  double t2862;
  double t2863;
  double t2854;
  double t2856;
  double t2857;
  double t2868;
  double t2869;
  double t3632;
  double t6498;
  double t6508;
  double t6512;
  double t2851;
  double t2864;
  double t6569;
  double t6576;
  double t6588;
  double t6596;
  double t8730;
  double t8749;
  double t8769;
  double t8776;
  double t8777;
  double t8780;
  double t8803;
  double t8804;
  double t8805;
  double t8846;
  double t8847;
  double t8809;
  double t8810;
  double t8811;
  double t8812;
  double t8826;
  double t8827;
  double t8828;
  double t8829;
  double t8830;
  double t8831;
  double t8832;
  double t8833;
  double t8834;
  double t8835;
  double t8806;
  double t8807;
  double t8808;
  double t8813;
  double t8818;
  double t8836;
  double t8837;
  double t8838;
  double t8839;
  double t8840;
  double t8841;
  double t8842;
  double t8843;
  double t2821;
  double t2832;
  double t2846;
  double t2858;
  double t3633;
  double t6476;
  double t6513;
  double t6521;
  double t6525;
  double t2828;
  double t2839;
  double t2852;
  double t2865;
  double t8844;
  double t6526;
  double t8891;
  double t8893;
  double t8857;
  double t8892;
  double t8894;
  double t8863;
  double t8871;
  double t8872;
  double t8873;
  double t8874;
  double t8875;
  double t8866;
  double t8867;
  double t8848;
  double t8849;
  double t8850;
  double t8851;
  double t8741;
  double t8781;
  double t8784;
  double t8885;
  double t8886;
  double t8887;
  double t8888;
  double t8898;
  double t8899;
  double t8900;
  double t8923;
  double t8924;
  double t8876;
  double t8877;
  double t8878;
  double t8879;
  double t8916;
  double t8918;
  double t8802;
  double t8890;
  double t8922;
  double t8939;
  double t8940;
  double t8941;
  double t8880;
  double t8914;
  double t8937;
  double t8785;
  double t8889;
  double t8921;
  double t8942;
  double t8943;
  double t8944;
  double t8881;
  double t8915;
  double t8938;
  t2831 = Sin(var2[2]);
  t2820 = Cos(var2[2]);
  t2840 = Cos(var2[3]);
  t2842 = Sin(var2[3]);
  t2853 = Cos(var2[4]);
  t2855 = Sin(var2[4]);
  t2848 = t2820*t2840;
  t2849 = -1.*t2831*t2842;
  t2850 = t2848 + t2849;
  t2841 = t2840*t2831;
  t2844 = t2820*t2842;
  t2845 = t2841 + t2844;
  t2861 = t2820*t2853;
  t2862 = -1.*t2831*t2855;
  t2863 = t2861 + t2862;
  t2854 = t2853*t2831;
  t2856 = t2820*t2855;
  t2857 = t2854 + t2856;
  t2868 = -1.*t2840*t2831;
  t2869 = -1.*t2820*t2842;
  t3632 = t2868 + t2869;
  t6498 = -1.*t2853*t2831;
  t6508 = -1.*t2820*t2855;
  t6512 = t6498 + t6508;
  t2851 = Power(t2850,2);
  t2864 = Power(t2863,2);
  t6569 = Power(t2840,2);
  t6576 = -0.5*t6569;
  t6588 = Power(t2842,2);
  t6596 = -0.5*t6588;
  t8730 = t6576 + t6596;
  t8749 = Power(t2853,2);
  t8769 = -0.5*t8749;
  t8776 = Power(t2855,2);
  t8777 = -0.5*t8776;
  t8780 = t8769 + t8777;
  t8803 = -1.*var1[4];
  t8804 = var3[3] + t8803;
  t8805 = -2.5*t8804*t3632;
  t8846 = -1.*var1[2];
  t8847 = var3[2] + t8846;
  t8809 = -1.*var1[0];
  t8810 = var3[0] + t8809;
  t8811 = 10.*t3632*t2850;
  t8812 = 10.*t2845*t2850;
  t8826 = -1.*var1[1];
  t8827 = var3[1] + t8826;
  t8828 = Power(t3632,2);
  t8829 = 5.*t8828;
  t8830 = 5.*t3632*t2845;
  t8831 = 5.*t2851;
  t8832 = -1.*t2820*t2840;
  t8833 = t2831*t2842;
  t8834 = t8832 + t8833;
  t8835 = 5.*t2850*t8834;
  t8806 = -1.*var1[3];
  t8807 = var3[4] + t8806;
  t8808 = -2.5*t8807*t6512;
  t8813 = 10.*t6512*t2863;
  t8818 = 10.*t2857*t2863;
  t8836 = Power(t6512,2);
  t8837 = 5.*t8836;
  t8838 = 5.*t6512*t2857;
  t8839 = 5.*t2864;
  t8840 = -1.*t2820*t2853;
  t8841 = t2831*t2855;
  t8842 = t8840 + t8841;
  t8843 = 5.*t2863*t8842;
  t2821 = Power(t2820,2);
  t2832 = Power(t2831,2);
  t2846 = Power(t2845,2);
  t2858 = Power(t2857,2);
  t3633 = -5.*t3632*t2850;
  t6476 = -5.*t2845*t2850;
  t6513 = -5.*t6512*t2863;
  t6521 = -5.*t2857*t2863;
  t6525 = t3633 + t6476 + t6513 + t6521;
  t2828 = -25.*t2821;
  t2839 = -25.*t2832;
  t2852 = -5.*t2851;
  t2865 = -5.*t2864;
  t8844 = t8829 + t8830 + t8831 + t8835 + t8837 + t8838 + t8839 + t8843;
  t6526 = -10.*t2820;
  t8891 = -2.5*t8804*t8834;
  t8893 = 10.*t3632*t8834;
  t8857 = t8829 + t8830 + t8831 + t8835;
  t8892 = -2.5*t8807*t8842;
  t8894 = 10.*t6512*t8842;
  t8863 = t8837 + t8838 + t8839 + t8843;
  t8871 = 5.*t3632*t2850;
  t8872 = 5.*t2845*t2850;
  t8873 = 5.*t6512*t2863;
  t8874 = 5.*t2857*t2863;
  t8875 = t8871 + t8872 + t8873 + t8874;
  t8866 = 25.*t2821;
  t8867 = 25.*t2832;
  t8848 = -10.*t2831;
  t8849 = 5.*t3632*t8730;
  t8850 = 5.*t6512*t8780;
  t8851 = t8848 + t8849 + t8850;
  t8741 = -5.*t2850*t8730;
  t8781 = -5.*t2863*t8780;
  t8784 = t6526 + t8741 + t8781;
  t8885 = 10.*t2831;
  t8886 = -5.*t3632*t8730;
  t8887 = -5.*t6512*t8780;
  t8888 = t8885 + t8886 + t8887;
  t8898 = 5.*t8834*t8730;
  t8899 = 5.*t8842*t8780;
  t8900 = t6526 + t8898 + t8899;
  t8923 = -1.*var4[0]*t2845;
  t8924 = -1.*var4[1]*t2850;
  t8876 = 10.*t2820;
  t8877 = 5.*t2850*t8730;
  t8878 = 5.*t2863*t8780;
  t8879 = t8876 + t8877 + t8878;
  t8916 = Power(t8730,2);
  t8918 = Power(t8780,2);
  t8802 = 2.5*t2850;
  t8890 = 2.5*t3632;
  t8922 = 2.5*t8730;
  t8939 = -2.5*t8810*t3632;
  t8940 = -2.5*t8827*t8834;
  t8941 = t8939 + t8923 + t8924 + t8940;
  t8880 = -2.5*t2850;
  t8914 = -2.5*t3632;
  t8937 = -2.5*t8730;
  t8785 = 2.5*t2863;
  t8889 = 2.5*t6512;
  t8921 = 2.5*t8780;
  t8942 = -2.5*t8810*t6512;
  t8943 = -2.5*t8827*t8842;
  t8944 = t8942 + t8943;
  t8881 = -2.5*t2863;
  t8915 = -2.5*t6512;
  t8938 = -2.5*t8780;
  p_output1[0]=t2828 + t2839 - 5.*t2846 + t2852 - 5.*t2858 + t2865;
  p_output1[1]=t6525;
  p_output1[2]=t8784;
  p_output1[3]=t8785;
  p_output1[4]=t8802;
  p_output1[5]=t8805 + t8808 + t8810*(t8811 + t8812 + t8813 + t8818) + t8827*t8844 + t8847*t8851;
  p_output1[6]=t8805 + t8810*(t8811 + t8812) + 5.*t3632*t8730*t8847 + t8827*t8857;
  p_output1[7]=t8808 + t8810*(t8813 + t8818) + 5.*t6512*t8780*t8847 + t8827*t8863;
  p_output1[8]=5.*t2846 + 5.*t2858 + t8831 + t8839 + t8866 + t8867;
  p_output1[9]=t8875;
  p_output1[10]=t8879;
  p_output1[11]=t8880;
  p_output1[12]=t8881;
  p_output1[13]=-1.;
  p_output1[14]=t6525;
  p_output1[15]=t2828 + t2839 + t2852 + t2865 - 5.*t8828 - 5.*t8836;
  p_output1[16]=t8888;
  p_output1[17]=t8889;
  p_output1[18]=t8890;
  p_output1[19]=t8810*t8844 + t8891 + t8892 + t8827*(t8811 + t8813 + t8893 + t8894) + t8847*t8900;
  p_output1[20]=5.*t8730*t8834*t8847 + t8810*t8857 + t8891 + t8827*(t8811 + t8893);
  p_output1[21]=5.*t8780*t8842*t8847 + t8810*t8863 + t8892 + t8827*(t8813 + t8894);
  p_output1[22]=t8875;
  p_output1[23]=t8829 + t8831 + t8837 + t8839 + t8866 + t8867;
  p_output1[24]=t8851;
  p_output1[25]=t8914;
  p_output1[26]=t8915;
  p_output1[27]=-1.;
  p_output1[28]=t8784;
  p_output1[29]=t8888;
  p_output1[30]=-4.47 - 5.*t8916 - 5.*t8918;
  p_output1[31]=t8921;
  p_output1[32]=t8922;
  p_output1[33]=t8810*t8851 + t8827*t8900 + t8923 + t8924;
  p_output1[34]=5.*t3632*t8730*t8810 + 5.*t8730*t8827*t8834 + t8923 + t8924;
  p_output1[35]=5.*t6512*t8780*t8810 + 5.*t8780*t8827*t8842;
  p_output1[36]=t8879;
  p_output1[37]=t8851;
  p_output1[38]=4.47 + 5.*t8916 + 5.*t8918;
  p_output1[39]=t8937;
  p_output1[40]=t8938;
  p_output1[41]=t2850;
  p_output1[42]=t3632;
  p_output1[43]=t8802;
  p_output1[44]=t8890;
  p_output1[45]=t8922;
  p_output1[46]=-1.25;
  p_output1[47]=t8941;
  p_output1[48]=t8941;
  p_output1[49]=t8880;
  p_output1[50]=t8914;
  p_output1[51]=t8937;
  p_output1[52]=1.25;
  p_output1[53]=t2850;
  p_output1[54]=t3632;
  p_output1[55]=t8785;
  p_output1[56]=t8889;
  p_output1[57]=t8921;
  p_output1[58]=-1.25;
  p_output1[59]=t8944;
  p_output1[60]=t8944;
  p_output1[61]=t8881;
  p_output1[62]=t8915;
  p_output1[63]=t8938;
  p_output1[64]=1.25;
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

  double *var1,*var2,*var3,*var4;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 4)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Four input(s) required (var1,var2,var3,var4).");
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
  mrows = mxGetM(prhs[3]);
  ncols = mxGetN(prhs[3]);
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
    ( !(mrows == 2 && ncols == 1) && 
      !(mrows == 1 && ncols == 2))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var4 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
  var3 = mxGetPr(prhs[2]);
  var4 = mxGetPr(prhs[3]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 65, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2,var3,var4);


}

#else // MATLAB_MEX_FILE

#include "J_dxDiscreteMapSwingLegImpact.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_dxDiscreteMapSwingLegImpact_raw(double *p_output1, const double *var1,const double *var2,const double *var3,const double *var4)
{
  // Call Subroutines
  output1(p_output1, var1, var2, var3, var4);

}

}

#endif // MATLAB_MEX_FILE
