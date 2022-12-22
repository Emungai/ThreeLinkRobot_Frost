/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:41 GMT-04:00
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
  double t7738;
  double t7740;
  double t7737;
  double t7739;
  double t7741;
  double t7742;
  double t7743;
  double t7751;
  double t7752;
  double t7753;
  double t7745;
  double t7746;
  double t7747;
  double t7748;
  double t7749;
  double t7759;
  double t7760;
  double t7761;
  double t7755;
  double t7756;
  double t7757;
  double t7763;
  double t7764;
  double t7765;
  double t7744;
  double t7790;
  double t7791;
  double t7792;
  double t7793;
  double t7794;
  double t7754;
  double t7758;
  double t7769;
  double t7770;
  double t7771;
  double t7772;
  double t7773;
  double t7774;
  double t7775;
  double t7776;
  double t7777;
  double t7750;
  double t7796;
  double t7797;
  double t7798;
  double t7799;
  double t7800;
  double t7762;
  double t7766;
  double t7778;
  double t7779;
  double t7780;
  double t7781;
  double t7782;
  double t7783;
  double t7784;
  double t7785;
  double t7786;
  double t7787;
  double t7837;
  double t7839;
  double t7808;
  double t7838;
  double t7840;
  double t7814;
  double t7826;
  double t7827;
  double t7828;
  double t7829;
  double t7830;
  double t7817;
  double t7818;
  double t7819;
  double t7820;
  double t7789;
  double t7795;
  double t7801;
  double t7802;
  double t7844;
  double t7845;
  double t7846;
  double t7847;
  double t7831;
  double t7832;
  double t7833;
  double t7834;
  double t7879;
  double t7880;
  double t7881;
  double t7835;
  double t7861;
  double t7877;
  double t7882;
  double t7883;
  double t7884;
  double t7836;
  double t7862;
  double t7878;
  t7738 = Sin(var1[2]);
  t7740 = Cos(var1[2]);
  t7737 = Cos(var1[3]);
  t7739 = -1.*t7737*t7738;
  t7741 = Sin(var1[3]);
  t7742 = -1.*t7740*t7741;
  t7743 = t7739 + t7742;
  t7751 = t7740*t7737;
  t7752 = -1.*t7738*t7741;
  t7753 = t7751 + t7752;
  t7745 = Cos(var1[4]);
  t7746 = -1.*t7745*t7738;
  t7747 = Sin(var1[4]);
  t7748 = -1.*t7740*t7747;
  t7749 = t7746 + t7748;
  t7759 = t7740*t7745;
  t7760 = -1.*t7738*t7747;
  t7761 = t7759 + t7760;
  t7755 = t7737*t7738;
  t7756 = t7740*t7741;
  t7757 = t7755 + t7756;
  t7763 = t7745*t7738;
  t7764 = t7740*t7747;
  t7765 = t7763 + t7764;
  t7744 = 2.5*var2[3]*t7743;
  t7790 = Power(t7737,2);
  t7791 = -0.5*t7790;
  t7792 = Power(t7741,2);
  t7793 = -0.5*t7792;
  t7794 = t7791 + t7793;
  t7754 = -10.*t7743*t7753;
  t7758 = -10.*t7757*t7753;
  t7769 = Power(t7743,2);
  t7770 = -5.*t7769;
  t7771 = -5.*t7743*t7757;
  t7772 = Power(t7753,2);
  t7773 = -5.*t7772;
  t7774 = -1.*t7740*t7737;
  t7775 = t7738*t7741;
  t7776 = t7774 + t7775;
  t7777 = -5.*t7753*t7776;
  t7750 = 2.5*var2[4]*t7749;
  t7796 = Power(t7745,2);
  t7797 = -0.5*t7796;
  t7798 = Power(t7747,2);
  t7799 = -0.5*t7798;
  t7800 = t7797 + t7799;
  t7762 = -10.*t7749*t7761;
  t7766 = -10.*t7765*t7761;
  t7778 = Power(t7749,2);
  t7779 = -5.*t7778;
  t7780 = -5.*t7749*t7765;
  t7781 = Power(t7761,2);
  t7782 = -5.*t7781;
  t7783 = -1.*t7740*t7745;
  t7784 = t7738*t7747;
  t7785 = t7783 + t7784;
  t7786 = -5.*t7761*t7785;
  t7787 = t7770 + t7771 + t7773 + t7777 + t7779 + t7780 + t7782 + t7786;
  t7837 = 2.5*var2[3]*t7776;
  t7839 = -10.*t7743*t7776;
  t7808 = t7770 + t7771 + t7773 + t7777;
  t7838 = 2.5*var2[4]*t7785;
  t7840 = -10.*t7749*t7785;
  t7814 = t7779 + t7780 + t7782 + t7786;
  t7826 = -5.*t7743*t7753;
  t7827 = -5.*t7757*t7753;
  t7828 = -5.*t7749*t7761;
  t7829 = -5.*t7765*t7761;
  t7830 = t7826 + t7827 + t7828 + t7829;
  t7817 = Power(t7740,2);
  t7818 = -25.*t7817;
  t7819 = Power(t7738,2);
  t7820 = -25.*t7819;
  t7789 = 10.*t7738;
  t7795 = -5.*t7743*t7794;
  t7801 = -5.*t7749*t7800;
  t7802 = t7789 + t7795 + t7801;
  t7844 = 10.*t7740;
  t7845 = -5.*t7776*t7794;
  t7846 = -5.*t7785*t7800;
  t7847 = t7844 + t7845 + t7846;
  t7831 = -10.*t7740;
  t7832 = -5.*t7753*t7794;
  t7833 = -5.*t7761*t7800;
  t7834 = t7831 + t7832 + t7833;
  t7879 = 2.5*var2[0]*t7743;
  t7880 = 2.5*var2[1]*t7776;
  t7881 = t7879 + t7880;
  t7835 = 2.5*t7753;
  t7861 = 2.5*t7743;
  t7877 = 2.5*t7794;
  t7882 = 2.5*var2[0]*t7749;
  t7883 = 2.5*var2[1]*t7785;
  t7884 = t7882 + t7883;
  t7836 = 2.5*t7761;
  t7862 = 2.5*t7749;
  t7878 = 2.5*t7800;
  p_output1[0]=t7744 + t7750 + (t7754 + t7758 + t7762 + t7766)*var2[0] + t7787*var2[1] + t7802*var2[2];
  p_output1[1]=t7744 + (t7754 + t7758)*var2[0] + t7808*var2[1] - 5.*t7743*t7794*var2[2];
  p_output1[2]=t7750 + (t7762 + t7766)*var2[0] + t7814*var2[1] - 5.*t7749*t7800*var2[2];
  p_output1[3]=-5.*Power(t7757,2) - 5.*Power(t7765,2) + t7773 + t7782 + t7818 + t7820;
  p_output1[4]=t7830;
  p_output1[5]=t7834;
  p_output1[6]=t7835;
  p_output1[7]=t7836;
  p_output1[8]=t7837 + t7838 + t7787*var2[0] + (t7754 + t7762 + t7839 + t7840)*var2[1] + t7847*var2[2];
  p_output1[9]=t7837 + t7808*var2[0] + (t7754 + t7839)*var2[1] - 5.*t7776*t7794*var2[2];
  p_output1[10]=t7838 + t7814*var2[0] + (t7762 + t7840)*var2[1] - 5.*t7785*t7800*var2[2];
  p_output1[11]=t7830;
  p_output1[12]=t7770 + t7773 + t7779 + t7782 + t7818 + t7820;
  p_output1[13]=t7802;
  p_output1[14]=t7861;
  p_output1[15]=t7862;
  p_output1[16]=t7802*var2[0] + t7847*var2[1];
  p_output1[17]=-5.*t7743*t7794*var2[0] - 5.*t7776*t7794*var2[1];
  p_output1[18]=-5.*t7749*t7800*var2[0] - 5.*t7785*t7800*var2[1];
  p_output1[19]=t7834;
  p_output1[20]=t7802;
  p_output1[21]=-4.47 - 5.*Power(t7794,2) - 5.*Power(t7800,2);
  p_output1[22]=t7877;
  p_output1[23]=t7878;
  p_output1[24]=t7881;
  p_output1[25]=t7881;
  p_output1[26]=t7835;
  p_output1[27]=t7861;
  p_output1[28]=t7877;
  p_output1[29]=-1.25;
  p_output1[30]=t7884;
  p_output1[31]=t7884;
  p_output1[32]=t7836;
  p_output1[33]=t7862;
  p_output1[34]=t7878;
  p_output1[35]=-1.25;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 36, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_MmatDx_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_MmatDx_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
