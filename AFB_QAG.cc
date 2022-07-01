#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include "LHAPDF/LHAPDF.h"

using namespace LHAPDF;
using namespace std;

// SM input (v0 Mathematica)
// Constants
#define PI 3.14159265
#define GeVtofb 0.38937966e+12

// SM parameters
#define MZ 91.19
#define GammaZ 2.5
#define alphaEM 7.8125e-3
#define stheta2W 0.232035

// // SM input (v1 xFitter)
// // Constants
// #define PI 3.14159265
// #define GeVtofb 0.389379338e+12

// // SM parameters
// #define MZ 91.1876
// #define GammaZ 2.4952
// #define alphaEM 7.29735e-3
// #define stheta2W 0.23127

// // SM input (v2 MCFM)
// #define PI 3.14159265
// #define GeVtofb 0.38937966e+12
// 
// #define MZ 91.1876
// #define GammaZ 2.4952
// #define alphaEM 7.756247e-3
// #define stheta2W 0.2228972

// set collider energy and luminosity
#define energy 13000
#define lum 300

// define acceptance cuts
#define eta_cut 2.5
#define pT_cut 20

// define rapidity cuts
#define y_min 0.0
#define y_max 0.0 // y_max = 0 means no upper cut

// PDF set and grid
#define setname "CT14nnlo"
#define PDF_set 0

// PDFs section
// Initialize PDFs
void initLHAPDF();
const vector<LHAPDF::PDF*> pdf = mkPDFs(setname);

// Integration parameters
const int alloc_space = 1000;
gsl_integration_workspace * w = gsl_integration_workspace_alloc(alloc_space);
// Integration option
const int key = 2;
// Error of the integration
const double epsabs = 0;
const double epsrel = 1e-2;

double *propagators (double Minv)
{
    const double e = sqrt(4*PI*alphaEM);
    const double gsm = (e/(sqrt(stheta2W)*sqrt(1-stheta2W)))*sqrt(1+pow(stheta2W,2));
    
    const double smangle = atan(-stheta2W);
    
    // SM couplings
    static double *couplings_photon = new double[8];
    couplings_photon[0] = e*(2.0/3.0);
    couplings_photon[1] = 0;
    couplings_photon[2] = e*(-1.0/3.0);
    couplings_photon[3] = 0;
    couplings_photon[4] = e*(-1.0);
    couplings_photon[5] = 0;
    couplings_photon[6] = 0;
    couplings_photon[7] = 0;
    
    static double *couplings_Z = new double[8];
    couplings_Z[0] = (1.0/2.0)*gsm*(1.0/6.0)*(3*cos(smangle)+8*sin(smangle));
    couplings_Z[1] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
    couplings_Z[2] = (1.0/2.0)*gsm*(1.0/6.0)*(-3*cos(smangle)-4*sin(smangle));
    couplings_Z[3] = (1.0/2.0)*gsm*(-cos(smangle)/2.0);
    couplings_Z[4] = (1.0/2.0)*gsm*((-cos(smangle)/2.0)+(-2*sin(smangle)));
    couplings_Z[5] = (1.0/2.0)*gsm*(-cos(smangle)/2.0);
    couplings_Z[6] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
    couplings_Z[7] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
        
    // Even combination of couplings
    double even_photon_up = (pow(couplings_photon[0],2)+pow(couplings_photon[1],2))*(pow(couplings_photon[4],2)+pow(couplings_photon[5],2));
    double even_photon_down = (pow(couplings_photon[2],2)+pow(couplings_photon[3],2))*(pow(couplings_photon[4],2)+pow(couplings_photon[5],2));
    double even_interf_up = ((couplings_photon[0]*couplings_Z[0])+(couplings_photon[1]*couplings_Z[1]))*((couplings_photon[4]*couplings_Z[4])+(couplings_photon[5]*couplings_Z[5]));
    double even_interf_down = ((couplings_photon[2]*couplings_Z[2])+(couplings_photon[3]*couplings_Z[3]))*((couplings_photon[4]*couplings_Z[4])+(couplings_photon[5]*couplings_Z[5]));
    double even_Z_up = (pow(couplings_Z[0],2)+pow(couplings_Z[1],2))*(pow(couplings_Z[4],2)+pow(couplings_Z[5],2));
    double even_Z_down = (pow(couplings_Z[2],2)+pow(couplings_Z[3],2))*(pow(couplings_Z[4],2)+pow(couplings_Z[5],2));
    
    // Odd combination of couplings
    double odd_photon_up = 4*couplings_photon[0]*couplings_photon[1]*couplings_photon[4]*couplings_photon[5];
    double odd_photon_down = 4*couplings_photon[2]*couplings_photon[3]*couplings_photon[4]*couplings_photon[5];
    double odd_interf_up = (couplings_photon[0]*couplings_Z[1]+couplings_photon[1]*couplings_Z[0])*(couplings_photon[4]*couplings_Z[5]+couplings_photon[5]*couplings_Z[4]);
    double odd_interf_down = (couplings_photon[2]*couplings_Z[3]+couplings_photon[3]*couplings_Z[2])*(couplings_photon[4]*couplings_Z[5]+couplings_photon[5]*couplings_Z[4]);
    double odd_Z_up = 4*couplings_Z[0]*couplings_Z[1]*couplings_Z[4]*couplings_Z[5];
    double odd_Z_down = 4*couplings_Z[2]*couplings_Z[3]*couplings_Z[4]*couplings_Z[5];
    
    // Propagators squared and interference
    double photon_squared = 1.0/pow(Minv,4);
    double interference = 2.0*(-pow(Minv,2)*(pow(MZ,2)-pow(Minv,2)))/(pow(Minv,4)*((pow(pow(MZ,2)-pow(Minv,2),2))+pow(MZ,2)*pow(GammaZ,2)));
    double Z_squared = 1.0/(pow(pow(MZ,2)-pow(Minv,2),2)+pow(MZ,2)*pow(GammaZ,2));
    
    static double *propagators = new double[4];
    
    propagators[0] = (even_photon_up * photon_squared)+(even_interf_up * interference) + (even_Z_up * Z_squared);
    propagators[1] = (odd_photon_up * photon_squared)+(odd_interf_up * interference) + (odd_Z_up * Z_squared);
    propagators[2] = (even_photon_down * photon_squared)+(even_interf_down * interference) + (even_Z_down * Z_squared);
    propagators[3] = (odd_photon_down * photon_squared)+(odd_interf_down * interference) + (odd_Z_down * Z_squared);
        
    return propagators;
}

////UUBAR EVEN FORWARD Matrix element
double uubarEF_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
    
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1u = (pdf[PDF_set]->xfxQ(2, x1, Q))/x1;
    double f1c = (pdf[PDF_set]->xfxQ(4, x1, Q))/x1;
    double f2ubar = (pdf[PDF_set]->xfxQ(-2, x2, Q))/x2;
    double f2cbar = (pdf[PDF_set]->xfxQ(-4, x2, Q))/x2;

    // PDF combinations    
    double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(max(0., tanh(eta_cut-abs(y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0;
 
    double angular_integration_EF = (qqbar_cos_theta_max-qqbar_cos_theta_min)+(1.0/3.0)*(pow(qqbar_cos_theta_max,3)-pow(qqbar_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EF = dsigma*angular_integration_EF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // U-UBAR
    double uubarEF = uubar_PDF*dsigma_EF;
  
    double *propagator = propagators (Minv);
    
    return uubarEF * propagator[0];
}

////UUBAR EVEN FORWARD Integration in rapidity
double integration_uubarEF_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &uubarEF_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error);
    
    return 2*result;
}

////UUBAR EVEN FORWARD Integration in invariant mass
double integration_uubarEF (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_uubarEF_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////UUBAR EVEN BACKWARD Matrix element
double uubarEB_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1u = (pdf[PDF_set]->xfxQ(2, x1, Q))/x1;
    double f1c = (pdf[PDF_set]->xfxQ(4, x1, Q))/x1;
    double f2ubar = (pdf[PDF_set]->xfxQ(-2, x2, Q))/x2;
    double f2cbar = (pdf[PDF_set]->xfxQ(-4, x2, Q))/x2;

    // PDF combinations    
    double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
    
    // Angular integration limits
    double qbarq_cos_theta_max = 0;
    double qbarq_cos_theta_min = max(min(0., -tanh(eta_cut-abs(y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
 
    double angular_integration_EB = (qbarq_cos_theta_max-qbarq_cos_theta_min)+(1.0/3.0)*(pow(qbarq_cos_theta_max,3)-pow(qbarq_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EB = dsigma*angular_integration_EB;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // U-UBAR
    double uubarEB = uubar_PDF*dsigma_EB;
  
    double *propagator = propagators (Minv);
    return uubarEB * propagator[0];
}

////UUBAR EVEN BACKWARD Integration in rapidity
double integration_uubarEB_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &uubarEB_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////UUBAR EVEN BACKWARD Integration in invariant mass
double integration_uubarEB (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_uubarEB_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////UUBAR ODD FORWARD Matrix element
double uubarOF_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1u = (pdf[PDF_set]->xfxQ(2, x1, Q))/x1;
    double f1c = (pdf[PDF_set]->xfxQ(4, x1, Q))/x1;
    double f2ubar = (pdf[PDF_set]->xfxQ(-2, x2, Q))/x2;
    double f2cbar = (pdf[PDF_set]->xfxQ(-4, x2, Q))/x2;

    // PDF combinations    
    double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(max(0., tanh(eta_cut-abs(y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0;  
 
    double angular_integration_OF = pow(qqbar_cos_theta_max,2) - pow(qqbar_cos_theta_min,2);
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_OF = dsigma*angular_integration_OF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // U-UBAR
    double uubarOF = uubar_PDF*dsigma_OF;
  
    double *propagator = propagators (Minv);
    return uubarOF * propagator[1];
}

////UUBAR ODD FORWARD Integration in rapidity
double integration_uubarOF_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &uubarOF_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////UUBAR ODD FORWARD Integration in invariant mass
double integration_uubarOF (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_uubarOF_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////UUBAR ODD BACKWARD Matrix element
double uubarOB_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1u = (pdf[PDF_set]->xfxQ(2, x1, Q))/x1;
    double f1c = (pdf[PDF_set]->xfxQ(4, x1, Q))/x1;
    double f2ubar = (pdf[PDF_set]->xfxQ(-2, x2, Q))/x2;
    double f2cbar = (pdf[PDF_set]->xfxQ(-4, x2, Q))/x2;

    // PDF combinations    
    double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
    
    // Angular integration limits
    double qbarq_cos_theta_max = 0;
    double qbarq_cos_theta_min = max(min(0., -tanh(eta_cut-abs(y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
 
    double angular_integration_OB = pow(qbarq_cos_theta_max,2) - pow(qbarq_cos_theta_min,2);
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_OB = dsigma*angular_integration_OB;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // U-UBAR
    double uubarOB = uubar_PDF*dsigma_OB;
  
    double *propagator = propagators (Minv);
    return uubarOB * propagator[1];
}

////UUBAR ODD BACKWARD Integration in rapidity
double integration_uubarOB_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &uubarOB_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////UUBAR ODD BACKWARD Integration in invariant mass
double integration_uubarOB (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_uubarOB_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////UBARU EVEN FORWARD Matrix element
double ubaruEF_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1ubar = (pdf[PDF_set]->xfxQ(-2, x1, Q))/x1;
    double f1cbar = (pdf[PDF_set]->xfxQ(-4, x1, Q))/x1;
    double f2u = (pdf[PDF_set]->xfxQ(2, x2, Q))/x2;
    double f2c = (pdf[PDF_set]->xfxQ(4, x2, Q))/x2;

    // PDF combinations    
    double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
    
    // Angular integration limits
    double qbarq_cos_theta_max = 0;
    double qbarq_cos_theta_min = max(min(0., -tanh(eta_cut-abs(y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
 
    double angular_integration_EB = (qbarq_cos_theta_max-qbarq_cos_theta_min)+(1.0/3.0)*(pow(qbarq_cos_theta_max,3)-pow(qbarq_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EB = dsigma*angular_integration_EB;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // UBAR-U
    double ubaruEF = ubaru_PDF*dsigma_EB;
  
    double *propagator = propagators (Minv);
    return ubaruEF * propagator[0];
}

////UBARU EVEN FORWARD Integration in rapidity
double integration_ubaruEF_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &ubaruEF_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////UBARU EVEN FORWARD Integration in invariant mass
double integration_ubaruEF (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_ubaruEF_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////UBARU EVEN BACKWARD Matrix element
double ubaruEB_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1ubar = (pdf[PDF_set]->xfxQ(-2, x1, Q))/x1;
    double f1cbar = (pdf[PDF_set]->xfxQ(-4, x1, Q))/x1;
    double f2u = (pdf[PDF_set]->xfxQ(2, x2, Q))/x2;
    double f2c = (pdf[PDF_set]->xfxQ(4, x2, Q))/x2;

    // PDF combinations    
    double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(max(0., tanh(eta_cut-abs(y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0;  
 
    double angular_integration_EF = (qqbar_cos_theta_max-qqbar_cos_theta_min)+(1.0/3.0)*(pow(qqbar_cos_theta_max,3)-pow(qqbar_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EF = dsigma*angular_integration_EF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // UBAR-U
    double ubaruEB = ubaru_PDF*dsigma_EF;
  
    double *propagator = propagators (Minv);
    return ubaruEB * propagator[0];
}

////UBARU EVEN BACKWARD Integration in rapidity
double integration_ubaruEB_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &ubaruEB_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////UBARU EVEN BACKWARD Integration in invariant mass
double integration_ubaruEB (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_ubaruEB_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////UBARU ODD FORWARD Matrix element
double ubaruOF_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1ubar = (pdf[PDF_set]->xfxQ(-2, x1, Q))/x1;
    double f1cbar = (pdf[PDF_set]->xfxQ(-4, x1, Q))/x1;
    double f2u = (pdf[PDF_set]->xfxQ(2, x2, Q))/x2;
    double f2c = (pdf[PDF_set]->xfxQ(4, x2, Q))/x2;

    // PDF combinations    
    double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
    
    // Angular integration limits
    double qbarq_cos_theta_max = 0;
    double qbarq_cos_theta_min = max(min(0., -tanh(eta_cut-abs(y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
 
    double angular_integration_OB = pow(qbarq_cos_theta_max,2) - pow(qbarq_cos_theta_min,2);
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_OB = dsigma*angular_integration_OB;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // UBAR-U
    double ubaruOF = ubaru_PDF*dsigma_OB;
  
    double *propagator = propagators (Minv);
    return ubaruOF * propagator[1];
}
////UBARU ODD FORWARD Integration in rapidity
double integration_ubaruOF_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &ubaruOF_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////UBARU ODD FORWARD Integration in invariant mass
double integration_ubaruOF (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_ubaruOF_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////UBARU ODD BACKWARD Matrix element
double ubaruOB_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1ubar = (pdf[PDF_set]->xfxQ(-2, x1, Q))/x1;
    double f1cbar = (pdf[PDF_set]->xfxQ(-4, x1, Q))/x1;
    double f2u = (pdf[PDF_set]->xfxQ(2, x2, Q))/x2;
    double f2c = (pdf[PDF_set]->xfxQ(4, x2, Q))/x2;

    // PDF combinations    
    double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(max(0., tanh(eta_cut-abs(y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0; 
 
    double angular_integration_OF = pow(qqbar_cos_theta_max,2) - pow(qqbar_cos_theta_min,2);
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_OF = dsigma*angular_integration_OF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // UBAR-U
    double ubaruOB = ubaru_PDF*dsigma_OF;
  
    double *propagator = propagators (Minv);
    return ubaruOB * propagator[1];
}

////UBARU ODD BACKWARD Integration in rapidity
double integration_ubaruOB_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &ubaruOB_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////UBARU ODD BACKWARD Integration in invariant mass
double integration_ubaruOB (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_ubaruOB_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////DDBAR EVEN FORWARD Matrix element
double ddbarEF_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1d = (pdf[PDF_set]->xfxQ(1, x1, Q))/x1;
    double f1s = (pdf[PDF_set]->xfxQ(3, x1, Q))/x1;
    double f1b = (pdf[PDF_set]->xfxQ(5, x1, Q))/x1;
    double f2dbar = (pdf[PDF_set]->xfxQ(-1, x2, Q))/x2;
    double f2sbar = (pdf[PDF_set]->xfxQ(-3, x2, Q))/x2;
    double f2bbar = (pdf[PDF_set]->xfxQ(-5, x2, Q))/x2;

    // PDF combinations    
    double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(max(0., tanh(eta_cut-abs(y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0;
 
    double angular_integration_EF = (qqbar_cos_theta_max-qqbar_cos_theta_min)+(1.0/3.0)*(pow(qqbar_cos_theta_max,3)-pow(qqbar_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EF = dsigma*angular_integration_EF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // D-DBAR
    double ddbarEF = ddbar_PDF*dsigma_EF;
  
    double *propagator = propagators (Minv);
    return ddbarEF * propagator[2];
}

////DDBAR EVEN FORWARD Integration in rapidity
double integration_ddbarEF_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &ddbarEF_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////DDBAR EVEN FORWARD Integration in invariant mass
double integration_ddbarEF (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_ddbarEF_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////DDBAR EVEN BACKWARD Matrix element
double ddbarEB_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1d = (pdf[PDF_set]->xfxQ(1, x1, Q))/x1;
    double f1s = (pdf[PDF_set]->xfxQ(3, x1, Q))/x1;
    double f1b = (pdf[PDF_set]->xfxQ(5, x1, Q))/x1;
    double f2dbar = (pdf[PDF_set]->xfxQ(-1, x2, Q))/x2;
    double f2sbar = (pdf[PDF_set]->xfxQ(-3, x2, Q))/x2;
    double f2bbar = (pdf[PDF_set]->xfxQ(-5, x2, Q))/x2;

    // PDF combinations    
    double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
    
    // Angular integration limits
    double qbarq_cos_theta_max = 0;
    double qbarq_cos_theta_min = max(min(0., -tanh(eta_cut-abs(y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
 
    double angular_integration_EB = (qbarq_cos_theta_max-qbarq_cos_theta_min)+(1.0/3.0)*(pow(qbarq_cos_theta_max,3)-pow(qbarq_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EB = dsigma*angular_integration_EB;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // D-DBAR
    double ddbarEB = ddbar_PDF*dsigma_EB;
  
    double *propagator = propagators (Minv);
    return ddbarEB * propagator[2];
}

////DDBAR EVEN BACKWARD Integration in rapidity
double integration_ddbarEB_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &ddbarEB_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////DDBAR EVEN BACKWARD Integration in invariant mass
double integration_ddbarEB (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_ddbarEB_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////DDBAR ODD FORWARD Matrix element
double ddbarOF_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1d = (pdf[PDF_set]->xfxQ(1, x1, Q))/x1;
    double f1s = (pdf[PDF_set]->xfxQ(3, x1, Q))/x1;
    double f1b = (pdf[PDF_set]->xfxQ(5, x1, Q))/x1;
    double f2dbar = (pdf[PDF_set]->xfxQ(-1, x2, Q))/x2;
    double f2sbar = (pdf[PDF_set]->xfxQ(-3, x2, Q))/x2;
    double f2bbar = (pdf[PDF_set]->xfxQ(-5, x2, Q))/x2;

    // PDF combinations    
    double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(max(0., tanh(eta_cut-abs(y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0;  
 
    double angular_integration_OF = pow(qqbar_cos_theta_max,2) - pow(qqbar_cos_theta_min,2);
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_OF = dsigma*angular_integration_OF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // D-DBAR
    double ddbarOF = ddbar_PDF*dsigma_OF;
  
    double *propagator = propagators (Minv);
    return ddbarOF * propagator[3];
}

////DDBAR ODD FORWARD Integration in rapidity
double integration_ddbarOF_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &ddbarOF_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////DDBAR ODD FORWARD Integration in invariant mass
double integration_ddbarOF (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_ddbarOF_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////DDBAR ODD BACKWARD Matrix element
double ddbarOB_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1d = (pdf[PDF_set]->xfxQ(1, x1, Q))/x1;
    double f1s = (pdf[PDF_set]->xfxQ(3, x1, Q))/x1;
    double f1b = (pdf[PDF_set]->xfxQ(5, x1, Q))/x1;
    double f2dbar = (pdf[PDF_set]->xfxQ(-1, x2, Q))/x2;
    double f2sbar = (pdf[PDF_set]->xfxQ(-3, x2, Q))/x2;
    double f2bbar = (pdf[PDF_set]->xfxQ(-5, x2, Q))/x2;

    // PDF combinations    
    double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
    
    // Angular integration limits
    double qbarq_cos_theta_max = 0;
    double qbarq_cos_theta_min = max(min(0., -tanh(eta_cut-abs(y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
 
    double angular_integration_OB = pow(qbarq_cos_theta_max,2) - pow(qbarq_cos_theta_min,2);
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_OB = dsigma*angular_integration_OB;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // D-DBAR
    double ddbarOB = ddbar_PDF*dsigma_OB;
  
    double *propagator = propagators (Minv);
    return ddbarOB * propagator[3];
}

////DDBAR ODD BACKWARD Integration in rapidity
double integration_ddbarOB_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &ddbarOB_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////DDBAR ODD BACKWARD Integration in invariant mass
double integration_ddbarOB (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_ddbarOB_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////DBARD EVEN FORWARD Matrix element
double dbardEF_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1dbar = (pdf[PDF_set]->xfxQ(-1, x1, Q))/x1;
    double f1sbar = (pdf[PDF_set]->xfxQ(-3, x1, Q))/x1;
    double f1bbar = (pdf[PDF_set]->xfxQ(-5, x1, Q))/x1;
    double f2d = (pdf[PDF_set]->xfxQ(1, x2, Q))/x2;
    double f2s = (pdf[PDF_set]->xfxQ(3, x2, Q))/x2;
    double f2b = (pdf[PDF_set]->xfxQ(5, x2, Q))/x2;

    // PDF combinations    
    double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
    
    // Angular integration limits
    double qbarq_cos_theta_max = 0;
    double qbarq_cos_theta_min = max(min(0., -tanh(eta_cut-abs(y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
 
    double angular_integration_EB = (qbarq_cos_theta_max-qbarq_cos_theta_min)+(1.0/3.0)*(pow(qbarq_cos_theta_max,3)-pow(qbarq_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EB = dsigma*angular_integration_EB;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // DBAR-D
    double dbardEF = dbard_PDF*dsigma_EB;
  
    double *propagator = propagators (Minv);
    return dbardEF * propagator[2];
}

////DBARD EVEN FORWARD Integration in rapidity
double integration_dbardEF_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &dbardEF_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////DBARD EVEN FORWARD Integration in invariant mass
double integration_dbardEF (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_dbardEF_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////DBARD EVEN BACKWARD Matrix element
double dbardEB_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1dbar = (pdf[PDF_set]->xfxQ(-1, x1, Q))/x1;
    double f1sbar = (pdf[PDF_set]->xfxQ(-3, x1, Q))/x1;
    double f1bbar = (pdf[PDF_set]->xfxQ(-5, x1, Q))/x1;
    double f2d = (pdf[PDF_set]->xfxQ(1, x2, Q))/x2;
    double f2s = (pdf[PDF_set]->xfxQ(3, x2, Q))/x2;
    double f2b = (pdf[PDF_set]->xfxQ(5, x2, Q))/x2;

    // PDF combinations    
    double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(max(0., tanh(eta_cut-abs(y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0;  
 
    double angular_integration_EF = (qqbar_cos_theta_max-qqbar_cos_theta_min)+(1.0/3.0)*(pow(qqbar_cos_theta_max,3)-pow(qqbar_cos_theta_min,3));
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_EF = dsigma*angular_integration_EF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // DBAR-D
    double dbardEB = dbard_PDF*dsigma_EF;
  
    double *propagator = propagators (Minv);
    return dbardEB * propagator[2];
}

////DBARD EVEN BACKWARD Integration in rapidity
double integration_dbardEB_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &dbardEB_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////DBARD EVEN BACKWARD Integration in invariant mass
double integration_dbardEB (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_dbardEB_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////DBARD ODD FORWARD Matrix element
double dbardOF_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1dbar = (pdf[PDF_set]->xfxQ(-1, x1, Q))/x1;
    double f1sbar = (pdf[PDF_set]->xfxQ(-3, x1, Q))/x1;
    double f1bbar = (pdf[PDF_set]->xfxQ(-5, x1, Q))/x1;
    double f2d = (pdf[PDF_set]->xfxQ(1, x2, Q))/x2;
    double f2s = (pdf[PDF_set]->xfxQ(3, x2, Q))/x2;
    double f2b = (pdf[PDF_set]->xfxQ(5, x2, Q))/x2;

    // PDF combinations    
    double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
    
    // Angular integration limits
    double qbarq_cos_theta_max = 0;
    double qbarq_cos_theta_min = max(min(0., -tanh(eta_cut-abs(y))),-sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
 
    double angular_integration_OB = pow(qbarq_cos_theta_max,2) - pow(qbarq_cos_theta_min,2);
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_OB = dsigma*angular_integration_OB;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // DBAR-D
    double dbardOF = dbard_PDF*dsigma_OB;
  
    double *propagator = propagators (Minv);
    return dbardOF * propagator[3];
}

////DBARD ODD FORWARD Integration in rapidity
double integration_dbardOF_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &dbardOF_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////DBARD ODD FORWARD Integration in invariant mass
double integration_dbardOF (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_dbardOF_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}

////DBARD ODD BACKWARD Matrix element
double dbardOB_funct (double yreduced, void * params) {
    
    // Pass the invariant mass as parameter
    double Minv = *(double *) params;
        
    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    // Partons PDFs
    double f1dbar = (pdf[PDF_set]->xfxQ(-1, x1, Q))/x1;
    double f1sbar = (pdf[PDF_set]->xfxQ(-3, x1, Q))/x1;
    double f1bbar = (pdf[PDF_set]->xfxQ(-5, x1, Q))/x1;
    double f2d = (pdf[PDF_set]->xfxQ(1, x2, Q))/x2;
    double f2s = (pdf[PDF_set]->xfxQ(3, x2, Q))/x2;
    double f2b = (pdf[PDF_set]->xfxQ(5, x2, Q))/x2;

    // PDF combinations    
    double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
    
    // Angular integration limits
    double qqbar_cos_theta_max = min(max(0., tanh(eta_cut-abs(y))),sqrt(1-4*(pow(pT_cut,2)/pow(Minv,2))));
    double qqbar_cos_theta_min = 0; 
 
    double angular_integration_OF = pow(qqbar_cos_theta_max,2) - pow(qqbar_cos_theta_min,2);
    
    // Combination with angular integration (Forward - Backward for q-qbar)
    double dsigma_OF = dsigma*angular_integration_OF;

    // Covolution with PDFs (flipping direction for q-qbar & qbar-q)
    // DBAR-D
    double dbardOB = dbard_PDF*dsigma_OF;
  
    double *propagator = propagators (Minv);
    return dbardOB * propagator[3];
}

////DBARD ODD BACKWARD Integration in rapidity
double integration_dbardOB_y (double Minv, void * params) {
    
    double result, error;

    gsl_function F;
    F.function = &dbardOB_funct;
    F.params = &Minv;
    
    double inf = y_min / log(energy/Minv);
    double sup;
    if (y_max == 0.0) {
        sup = 1;
    } else {
        sup = y_max / log(energy/Minv);
    }
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return 2*result;
}

////DBARD ODD BACKWARD Integration in invariant mass
double integration_dbardOB (double Minv_inf, double Minv_sup) {
    
    double result, error;

    gsl_function F;
    F.function = &integration_dbardOB_y;
    
    double inf = Minv_inf;
    double sup = Minv_sup;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error); 

    return result;
}


double *observables (double Minv_inf, double Minv_sup)
{
    double uubarEF = integration_uubarEF (Minv_inf, Minv_sup);
    double uubarEB = integration_uubarEB (Minv_inf, Minv_sup);
    double uubarOF = integration_uubarOF (Minv_inf, Minv_sup);
    double uubarOB = integration_uubarOB (Minv_inf, Minv_sup);
    
    double ubaruEF = integration_ubaruEF (Minv_inf, Minv_sup);
    double ubaruEB = integration_ubaruEB (Minv_inf, Minv_sup);
    double ubaruOF = integration_ubaruOF (Minv_inf, Minv_sup);
    double ubaruOB = integration_ubaruOB (Minv_inf, Minv_sup);
    
    double ddbarEF = integration_ddbarEF (Minv_inf, Minv_sup);
    double ddbarEB = integration_ddbarEB (Minv_inf, Minv_sup);
    double ddbarOF = integration_ddbarOF (Minv_inf, Minv_sup);
    double ddbarOB = integration_ddbarOB (Minv_inf, Minv_sup);
        
    double dbardEF = integration_dbardEF (Minv_inf, Minv_sup); 
    double dbardEB = integration_dbardEB (Minv_inf, Minv_sup);
    double dbardOF = integration_dbardOF (Minv_inf, Minv_sup);
    double dbardOB = integration_dbardOB (Minv_inf, Minv_sup);
    
    // Reconstructed Forward and Backward
    double Forward = uubarEF+ubaruEF+ddbarEF+dbardEF+uubarOF+ubaruOF+ddbarOF+dbardOF;
    double Backward = uubarEB+ubaruEB+ddbarEB+dbardEB+uubarOB+ubaruOB+ddbarOB+dbardOB;

    // Cross section and statistical relative uncertainty
    double XS = Forward + Backward;
    double Epsilon_XS = 100.0 / sqrt(XS*lum);
    
    
    // Reconstructed AFB and statistical absolute uncertainty
    double AFB = (Forward - Backward) / (Forward + Backward);
    double Delta_AFB = sqrt((1.0-pow(AFB,2)) / (XS*lum));
    
    
    double *results = new double[4];
    
    results[0] = XS;
    results[1] = Epsilon_XS;
    results[2] = AFB;
    results[3] = Delta_AFB;
    
    return results;
}

int main(int argc, char** argv)
{
    // Value of the invariant mass is passed from command line
    if (argc <= 2)
    {
        printf ("Not enough entries, please enter the invariant mass\n");
        return 0;
    }    
    
    double Minv_inf, Minv_sup;
    Minv_inf = strtod (argv[1], NULL);
    Minv_sup = strtod (argv[2], NULL);
    
    // check on the rapidity cut
    if (y_min >= eta_cut) {
        printf("\nThe chosen lower rapidity cut is higher than pseudorapidity cut.\n\n");
        return 0;
    }
    if (y_min / log(energy/Minv_sup) > 1) {
        printf("\nThe chosen lower rapidity cut is too high in this invariant mass range.\n\n");
        return 0;
    }
        
    double *results = observables (Minv_inf, Minv_sup);
    
    gsl_integration_workspace_free (w);
    
    double XS_result = results[0];
    double Epsilon_XS_result = results[1];
    double AFB_result = results[2];
    double Delta_AFB_result = results[3];
    
    printf ("XS [fb]          = % .18f\n", XS_result);
    printf ("\u03B5XS [%]          = % .18f\n", Epsilon_XS_result);
    printf ("AFB*          = % .18f\n", AFB_result);
    printf ("\u0394AFB*          = % .18f\n", Delta_AFB_result);

    return 0;
}
