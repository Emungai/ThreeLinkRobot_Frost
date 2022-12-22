/*
 * Automatically Generated from Mathematica.
 * Wed 2 Aug 2017 17:44:49 GMT-04:00
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
  double t8063;
  double t8071;
  double t8064;
  double t8080;
  double t8093;
  double t8094;
  double t8095;
  double t8090;
  double t8091;
  double t8092;
  double t8065;
  double t8081;
  double t8082;
  double t8097;
  double t8098;
  double t8099;
  double t8084;
  double t8086;
  double t8106;
  double t8107;
  double t8108;
  double t8103;
  double t8104;
  double t8105;
  double t8085;
  double t8087;
  double t8088;
  double t8110;
  double t8111;
  double t8112;
  double t8083;
  double t8133;
  double t8134;
  double t8135;
  double t8136;
  double t8137;
  double t8096;
  double t8100;
  double t8101;
  double t8102;
  double t8118;
  double t8119;
  double t8120;
  double t8121;
  double t8122;
  double t8123;
  double t8089;
  double t8139;
  double t8140;
  double t8141;
  double t8142;
  double t8143;
  double t8109;
  double t8113;
  double t8114;
  double t8115;
  double t8124;
  double t8125;
  double t8126;
  double t8127;
  double t8128;
  double t8129;
  double t8163;
  double t8164;
  double t8165;
  double t8166;
  double t8167;
  double t8169;
  double t8170;
  double t8171;
  double t8172;
  double t8173;
  double t8174;
  double t8175;
  double t8176;
  double t8177;
  double t8116;
  double t8191;
  double t8150;
  double t8194;
  double t8195;
  double t8192;
  double t8157;
  double t8196;
  double t8197;
  double t8178;
  double t8219;
  double t8220;
  double t8221;
  double t8132;
  double t8138;
  double t8144;
  double t8145;
  double t8200;
  double t8201;
  double t8202;
  double t8203;
  double t8183;
  double t8184;
  double t8185;
  double t8186;
  double t8146;
  double t8245;
  double t8246;
  double t8247;
  double t8248;
  double t8189;
  double t8227;
  double t8252;
  double t8253;
  double t8254;
  double t8255;
  double t8190;
  double t8228;
  t8063 = Cos(var1[2]);
  t8071 = Sin(var1[2]);
  t8064 = Cos(var1[3]);
  t8080 = Sin(var1[3]);
  t8093 = t8063*t8064;
  t8094 = -1.*t8071*t8080;
  t8095 = t8093 + t8094;
  t8090 = -1.*t8064*t8071;
  t8091 = -1.*t8063*t8080;
  t8092 = t8090 + t8091;
  t8065 = -1.*t8063*t8064;
  t8081 = t8071*t8080;
  t8082 = t8065 + t8081;
  t8097 = t8064*t8071;
  t8098 = t8063*t8080;
  t8099 = t8097 + t8098;
  t8084 = Cos(var1[4]);
  t8086 = Sin(var1[4]);
  t8106 = t8063*t8084;
  t8107 = -1.*t8071*t8086;
  t8108 = t8106 + t8107;
  t8103 = -1.*t8084*t8071;
  t8104 = -1.*t8063*t8086;
  t8105 = t8103 + t8104;
  t8085 = -1.*t8063*t8084;
  t8087 = t8071*t8086;
  t8088 = t8085 + t8087;
  t8110 = t8084*t8071;
  t8111 = t8063*t8086;
  t8112 = t8110 + t8111;
  t8083 = 1.25*var2[3]*t8082;
  t8133 = Power(t8064,2);
  t8134 = -0.5*t8133;
  t8135 = Power(t8080,2);
  t8136 = -0.5*t8135;
  t8137 = t8134 + t8136;
  t8096 = 15.*t8092*t8095;
  t8100 = 5.*t8099*t8095;
  t8101 = 15.*t8092*t8082;
  t8102 = 5.*t8099*t8082;
  t8118 = Power(t8092,2);
  t8119 = 10.*t8118;
  t8120 = 10.*t8092*t8099;
  t8121 = Power(t8095,2);
  t8122 = 10.*t8121;
  t8123 = 10.*t8095*t8082;
  t8089 = 1.25*var2[4]*t8088;
  t8139 = Power(t8084,2);
  t8140 = -0.5*t8139;
  t8141 = Power(t8086,2);
  t8142 = -0.5*t8141;
  t8143 = t8140 + t8142;
  t8109 = 15.*t8105*t8108;
  t8113 = 5.*t8112*t8108;
  t8114 = 15.*t8105*t8088;
  t8115 = 5.*t8112*t8088;
  t8124 = Power(t8105,2);
  t8125 = 10.*t8124;
  t8126 = 10.*t8105*t8112;
  t8127 = Power(t8108,2);
  t8128 = 10.*t8127;
  t8129 = 10.*t8108*t8088;
  t8163 = 10.*t8092*t8095;
  t8164 = 10.*t8099*t8095;
  t8165 = 10.*t8105*t8108;
  t8166 = 10.*t8112*t8108;
  t8167 = t8163 + t8164 + t8165 + t8166;
  t8169 = 5.*t8118;
  t8170 = 5.*t8092*t8099;
  t8171 = 5.*t8121;
  t8172 = 5.*t8095*t8082;
  t8173 = 5.*t8124;
  t8174 = 5.*t8105*t8112;
  t8175 = 5.*t8127;
  t8176 = 5.*t8108*t8088;
  t8177 = t8169 + t8170 + t8171 + t8172 + t8173 + t8174 + t8175 + t8176;
  t8116 = t8096 + t8100 + t8101 + t8102 + t8109 + t8113 + t8114 + t8115;
  t8191 = 1.25*var2[3]*t8099;
  t8150 = t8096 + t8100 + t8101 + t8102;
  t8194 = Power(t8082,2);
  t8195 = 10.*t8194;
  t8192 = 1.25*var2[4]*t8112;
  t8157 = t8109 + t8113 + t8114 + t8115;
  t8196 = Power(t8088,2);
  t8197 = 10.*t8196;
  t8178 = -0.5*var2[2]*t8177;
  t8219 = 10.*t8092*t8082;
  t8220 = 10.*t8105*t8088;
  t8221 = t8163 + t8219 + t8165 + t8220;
  t8132 = -10.*t8063;
  t8138 = 5.*t8082*t8137;
  t8144 = 5.*t8088*t8143;
  t8145 = t8132 + t8138 + t8144;
  t8200 = 10.*t8071;
  t8201 = 5.*t8099*t8137;
  t8202 = 5.*t8112*t8143;
  t8203 = t8200 + t8201 + t8202;
  t8183 = -10.*t8071;
  t8184 = 5.*t8092*t8137;
  t8185 = 5.*t8105*t8143;
  t8186 = t8183 + t8184 + t8185;
  t8146 = -0.5*var2[2]*t8145;
  t8245 = 1.25*var2[1]*t8099;
  t8246 = 1.25*var2[0]*t8082;
  t8247 = t8245 + t8246;
  t8248 = var2[2]*t8247;
  t8189 = 1.25*var2[2]*t8092;
  t8227 = 1.25*var2[2]*t8082;
  t8252 = 1.25*var2[1]*t8112;
  t8253 = 1.25*var2[0]*t8088;
  t8254 = t8252 + t8253;
  t8255 = var2[2]*t8254;
  t8190 = 1.25*var2[2]*t8105;
  t8228 = 1.25*var2[2]*t8088;
  p_output1[0]=(t8083 + t8089 + t8146 - 0.5*(t8119 + t8120 + t8122 + t8123 + t8125 + t8126 + t8128 + t8129)*var2[0] - 0.5*t8116*var2[1])*var2[2];
  p_output1[1]=var2[2]*(t8083 - 0.5*(t8119 + t8120 + t8122 + t8123)*var2[0] - 0.5*t8150*var2[1] - 2.5*t8082*t8137*var2[2]);
  p_output1[2]=var2[2]*(t8089 - 0.5*(t8125 + t8126 + t8128 + t8129)*var2[0] - 0.5*t8157*var2[1] - 2.5*t8088*t8143*var2[2]);
  p_output1[3]=-0.5*t8167*var2[2];
  p_output1[4]=t8178;
  p_output1[5]=-0.5*t8167*var2[0] - 0.5*t8177*var2[1] - 1.*t8186*var2[2] + 1.25*t8092*var2[3] + 1.25*t8105*var2[4];
  p_output1[6]=t8189;
  p_output1[7]=t8190;
  p_output1[8]=var2[2]*(t8191 + t8192 - 0.5*t8116*var2[0] - 0.5*(t8119 + t8120 + t8123 + t8125 + t8126 + t8129 + t8195 + t8197)*var2[1] - 0.5*t8203*var2[2]);
  p_output1[9]=var2[2]*(t8191 - 0.5*t8150*var2[0] - 0.5*(t8119 + t8120 + t8123 + t8195)*var2[1] - 2.5*t8099*t8137*var2[2]);
  p_output1[10]=var2[2]*(t8192 - 0.5*t8157*var2[0] - 0.5*(t8125 + t8126 + t8129 + t8197)*var2[1] - 2.5*t8112*t8143*var2[2]);
  p_output1[11]=t8178;
  p_output1[12]=-0.5*t8221*var2[2];
  p_output1[13]=t8083 + t8089 - 0.5*t8177*var2[0] - 0.5*t8221*var2[1] - 1.*t8145*var2[2];
  p_output1[14]=t8227;
  p_output1[15]=t8228;
  p_output1[16]=(-0.5*t8145*var2[0] - 0.5*t8203*var2[1])*var2[2];
  p_output1[17]=(-2.5*t8082*t8137*var2[0] - 2.5*t8099*t8137*var2[1])*var2[2];
  p_output1[18]=(-2.5*t8088*t8143*var2[0] - 2.5*t8112*t8143*var2[1])*var2[2];
  p_output1[19]=-0.5*t8186*var2[2];
  p_output1[20]=t8146;
  p_output1[21]=-0.5*t8186*var2[0] - 0.5*t8145*var2[1];
  p_output1[22]=t8248;
  p_output1[23]=t8248;
  p_output1[24]=t8189;
  p_output1[25]=t8227;
  p_output1[26]=1.25*t8092*var2[0] + 1.25*t8082*var2[1];
  p_output1[27]=t8255;
  p_output1[28]=t8255;
  p_output1[29]=t8190;
  p_output1[30]=t8228;
  p_output1[31]=1.25*t8105*var2[0] + 1.25*t8088*var2[1];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 32, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "J_Ce2_vec3_Three_link_walker.hh"

namespace Pattern[ThreeLink, Blank[system]]
{

void J_Ce2_vec3_Three_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
