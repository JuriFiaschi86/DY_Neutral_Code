#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <regex>
#include <chrono>
#include <gsl/gsl_integration.h>
#include "LHAPDF/LHAPDF.h"

using namespace LHAPDF;
using namespace std;

// Constants
#define PI 3.14159265
#define GeVtofb 0.38937966e+12

// set collider energy and luminosity
#define energy 13000

// define acceptance cuts
#define eta_cut 2.5
#define pT_cut 20

// #define eta_cut 100
// #define pT_cut 0

// define rapidity cuts
#define y_min 0.0
#define y_max 0.0 // y_max = 0 means no upper cut

const int PDF_error_flag = 0;

// PDF set and grid
// #define setname "CT18NNLO68cl"

// #define setname "CT18NNLO68cl_AFB_300"
// #define setname "CT18NNLO68cl_AW_300"
// #define setname "CT18NNLO68cl_AFBonAW_300"
// #define setname "CT18NNLO68cl_AFB_3000"
// #define setname "CT18NNLO68cl_AW_3000"
// #define setname "CT18NNLO68cl_AFBonAW_3000"

// #define setname "CT18NNLO68cl_AFB_300_hm"
// #define setname "CT18NNLO68cl_AW_300_hm"
// #define setname "CT18NNLO68cl_AFBonAW_300_hm"
// #define setname "CT18NNLO68cl_AFB_3000_hm"
// #define setname "CT18NNLO68cl_AW_3000_hm"
// #define setname "CT18NNLO68cl_AFBonAW_3000_hm"

// #define setname "CT18NNLO68cl_AWonAFB_300_hm"

// #define setname "MMHT2014nnlo68cl"

// #define setname "CT18NNLO68cl_AW_40_540_300"
// #define setname "CT18NNLO68cl_AW_500_1000_300"
// #define setname "CT18NNLO68cl_AW_1000_1500_300"
// #define setname "CT18NNLO68cl_AW_1500_2000_300"
// #define setname "CT18NNLO68cl_AW_40_540_3000"
// #define setname "CT18NNLO68cl_AW_500_1000_3000"
// #define setname "CT18NNLO68cl_AW_1000_1500_3000"
// #define setname "CT18NNLO68cl_AW_1500_2000_3000"

// #define setname "CT18NNLO68cl_AFB_250_750_300"
// #define setname "CT18NNLO68cl_AFB_1000_1500_300"
// #define setname "CT18NNLO68cl_AFB_1500_2000_300"
// #define setname "CT18NNLO68cl_AFB_250_750_3000"
// #define setname "CT18NNLO68cl_AFB_1000_1500_3000"
// #define setname "CT18NNLO68cl_AFB_1500_2000_3000"




// #define setname "MSHT20nnlo_as118"
// #define setname "MSHT20nnlo_as118_AFB_300"
// #define setname "MSHT20nnlo_as118_AFB_250_750_300"
// #define setname "MSHT20nnlo_as118_AFB_500_1000_300"
// #define setname "MSHT20nnlo_as118_AFB_1000_1500_300"
// #define setname "MSHT20nnlo_as118_AFB_1500_2000_300"
#define setname "MSHT20nnlo_as118_AFB_3000"
// #define setname "MSHT20nnlo_as118_AFB_250_750_3000"
// #define setname "MSHT20nnlo_as118_AFB_500_1000_3000"
// #define setname "MSHT20nnlo_as118_AFB_1000_1500_3000"
// #define setname "MSHT20nnlo_as118_AFB_1500_2000_3000"

// #define setname "HERAPDF20_NNLO_EIG"
// #define setname "HERAPDF20_NNLO_EIG_AFB_300"
// #define setname "HERAPDF20_NNLO_EIG_AFB_250_750_300"
// #define setname "HERAPDF20_NNLO_EIG_AFB_500_1000_300"
// #define setname "HERAPDF20_NNLO_EIG_AFB_1000_1500_300"
// #define setname "HERAPDF20_NNLO_EIG_AFB_1500_2000_300"
// #define setname "HERAPDF20_NNLO_EIG_AFB_3000"
// #define setname "HERAPDF20_NNLO_EIG_AFB_250_750_3000"
// #define setname "HERAPDF20_NNLO_EIG_AFB_500_1000_3000"
// #define setname "HERAPDF20_NNLO_EIG_AFB_1000_1500_3000"
// #define setname "HERAPDF20_NNLO_EIG_AFB_1500_2000_3000"


// #define setname "NNPDF31_nnlo_as_0118_hessian"
// #define setname "NNPDF31_nnlo_as_0118_hessian_AFB_300"
// #define setname "NNPDF31_nnlo_as_0118_hessian_AFB_250_750_300"
// #define setname "NNPDF31_nnlo_as_0118_hessian_AFB_500_1000_300"
// #define setname "NNPDF31_nnlo_as_0118_hessian_AFB_1000_1500_300"
// #define setname "NNPDF31_nnlo_as_0118_hessian_AFB_1500_2000_300"






const LHAPDF::PDFSet set(setname);
const size_t nmem = set.size()-1;
// const size_t nmem = 0;

const vector<LHAPDF::PDF*> pdf = set.mkPDFs();

// Invariant mass range and number of points
const double Minv_min = 2*pT_cut + 0.001;
const double Minv_max = energy;
const int points_Minv = 1200;
const double step = (Minv_max - Minv_min) / points_Minv;

// Integration parameters
const int alloc_space = points_Minv;
gsl_integration_workspace * w = gsl_integration_workspace_alloc(alloc_space);
// Integration option
const int key = 6;
// Error of the integration
const double epsabs = 0;
const double epsrel = 1e-2;


struct param_struct {
    double Minv;
    int PDF_set;
};


////UUBAR EVEN FORWARD Matrix element
double uubarEF_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
    
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
    
    return uubarEF;
}


////UUBAR EVEN FORWARD Integration in rapidity
double integration_uubarEF_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &uubarEF_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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


////UUBAR EVEN BACKWARD Matrix element
double uubarEB_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return uubarEB;
}

////UUBAR EVEN BACKWARD Integration in rapidity
double integration_uubarEB_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &uubarEB_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////UUBAR ODD FORWARD Matrix element
double uubarOF_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return uubarOF;
}

////UUBAR ODD FORWARD Integration in rapidity
double integration_uubarOF_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &uubarOF_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////UUBAR ODD BACKWARD Matrix element
double uubarOB_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return uubarOB;
}

////UUBAR ODD BACKWARD Integration in rapidity
double integration_uubarOB_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &uubarOB_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////UBARU EVEN FORWARD Matrix element
double ubaruEF_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return ubaruEF;
}

////UBARU EVEN FORWARD Integration in rapidity
double integration_ubaruEF_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &ubaruEF_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////UBARU EVEN BACKWARD Matrix element
double ubaruEB_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return ubaruEB;
}

////UBARU EVEN BACKWARD Integration in rapidity
double integration_ubaruEB_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &ubaruEB_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////UBARU ODD FORWARD Matrix element
double ubaruOF_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return ubaruOF;
}
////UBARU ODD FORWARD Integration in rapidity
double integration_ubaruOF_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &ubaruOF_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////UBARU ODD BACKWARD Matrix element
double ubaruOB_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return ubaruOB;
}

////UBARU ODD BACKWARD Integration in rapidity
double integration_ubaruOB_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &ubaruOB_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////DDBAR EVEN FORWARD Matrix element
double ddbarEF_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return ddbarEF;
}

////DDBAR EVEN FORWARD Integration in rapidity
double integration_ddbarEF_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &ddbarEF_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////DDBAR EVEN BACKWARD Matrix element
double ddbarEB_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return ddbarEB;
}

////DDBAR EVEN BACKWARD Integration in rapidity
double integration_ddbarEB_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &ddbarEB_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////DDBAR ODD FORWARD Matrix element
double ddbarOF_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return ddbarOF;
}

////DDBAR ODD FORWARD Integration in rapidity
double integration_ddbarOF_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &ddbarOF_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////DDBAR ODD BACKWARD Matrix element
double ddbarOB_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return ddbarOB;
}

////DDBAR ODD BACKWARD Integration in rapidity
double integration_ddbarOB_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &ddbarOB_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////DBARD EVEN FORWARD Matrix element
double dbardEF_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return dbardEF;
}

////DBARD EVEN FORWARD Integration in rapidity
double integration_dbardEF_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &dbardEF_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////DBARD EVEN BACKWARD Matrix element
double dbardEB_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return dbardEB;
}

////DBARD EVEN BACKWARD Integration in rapidity
double integration_dbardEB_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &dbardEB_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////DBARD ODD FORWARD Matrix element
double dbardOF_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return dbardOF;
}

////DBARD ODD FORWARD Integration in rapidity
double integration_dbardOF_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &dbardOF_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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

////DBARD ODD BACKWARD Matrix element
double dbardOB_funct (double yreduced, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double Minv = (p->Minv);
    int PDF_set = (p->PDF_set);
        
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
  
    return dbardOB;
}

////DBARD ODD BACKWARD Integration in rapidity
double integration_dbardOB_y (int PDF_set, double Minv) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &dbardOB_funct;
    F.params = &p;
    
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
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



int main(int argc, char** argv)
{
 
    // Time counter
    std::chrono::steady_clock::time_point begin, end;
    
    // Output file name
    char filename[80];
    strcpy(filename, "DY_");
    strcat(filename, setname);
    strcat(filename, "_eigen_v2.dat");
    
    ofstream out(filename);    
    
    // check on the rapidity cut
    if (y_min >= eta_cut) {
        printf("\nThe chosen lower rapidity cut is higher than pseudorapidity cut.\n\n");
        return 0;
    }
    if (y_min / log(energy/Minv_max) > 1) {
        printf("\nThe chosen lower rapidity cut is too high in this invariant mass range.\n\n");
        return 0;
    }
    
    vector<double> Minv_table;
    // Invariant mass table
    Minv_table.push_back(Minv_min);
    for (int i = 1; i < points_Minv; i++) {
        Minv_table.push_back(Minv_table[i-1]+step);
    }
    
    // Set strings
    stringstream eta_cut_string;
    eta_cut_string << fixed << setprecision(1) << eta_cut;
    string eta_cut_s = eta_cut_string.str();
    
    stringstream pT_cut_string;
    pT_cut_string << fixed << setprecision(0) << pT_cut;
    string pT_cut_s = pT_cut_string.str();
    
    stringstream energy_string;
    energy_string << fixed << setprecision(0) << energy;
    string energy_s = energy_string.str();
    
    string text;
    
    // UUBAR EVEN FORWARD
    cout << "Computing UUBAR EVEN FORWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double uubarEF;
    text = "listuubcutEF[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            uubarEF = integration_uubarEF_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << uubarEF;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UUBAR EVEN FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;

    
    // UUBAR EVEN BACKWARD
    cout << "Computing UUBAR EVEN BACKWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double uubarEB;
    text = "listuubcutEB[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            uubarEB = integration_uubarEB_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << uubarEB;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UUBAR EVEN BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // UUBAR ODD FORWARD
    cout << "Computing UUBAR ODD FORWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double uubarOF;
    text = "listuubcutOF[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            uubarOF = integration_uubarOF_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << uubarOF;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UUBAR ODD FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // UUBAR ODD BACKWARD
    cout << "Computing UUBAR ODD BACKWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double uubarOB;
    text = "listuubcutOB[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            uubarOB = integration_uubarOB_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << uubarOB;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UUBAR ODD BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // UBARU EVEN FORWARD
    cout << "Computing UBARU EVEN FORWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double ubaruEF;
    text = "listubucutEF[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            ubaruEF = integration_ubaruEF_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << ubaruEF;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UBARU EVEN FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // UBARU EVEN BACKWARD
    cout << "Computing UBARU EVEN BACKWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double ubaruEB;
    text = "listubucutEB[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            ubaruEB = integration_ubaruEB_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << ubaruEB;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UBARU EVEN BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // UBARU ODD FORWARD
    cout << "Computing UBARU ODD FORWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double ubaruOF;
    text = "listubucutOF[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            ubaruOF = integration_ubaruOF_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << ubaruOF;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UBARU ODD FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // UBARU ODD BACKWARD
    cout << "Computing UBARU ODD BACKWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double ubaruOB;
    text = "listubucutOB[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            ubaruOB = integration_ubaruOB_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << ubaruOB;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UBARU ODD BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // DDBAR EVEN FORWARD
    cout << "Computing DDBAR EVEN FORWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double ddbarEF;
    text = "listddbcutEF[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            ddbarEF = integration_ddbarEF_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << ddbarEF;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DDBAR EVEN FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // DDBAR EVEN BACKWARD
    cout << "Computing DDBAR EVEN BACKWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double ddbarEB;
    text = "listddbcutEB[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            ddbarEB = integration_ddbarEB_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << ddbarEB;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DDBAR EVEN BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // DDBAR ODD FORWARD
    cout << "Computing DDBAR ODD FORWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double ddbarOF;
    text = "listddbcutOF[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            ddbarOF = integration_ddbarOF_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << ddbarOF;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DDBAR ODD FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // DDBAR ODD BACKWARD
    cout << "Computing DDBAR ODD BACKWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double ddbarOB;
    text = "listddbcutOB[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            ddbarOB = integration_ddbarOB_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << ddbarOB;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DDBAR ODD BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // DBARD EVEN FORWARD
    cout << "Computing DBARD EVEN FORWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double dbardEF;
    text = "listdbdcutEF[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            dbardEF = integration_dbardEF_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << dbardEF;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DBARD EVEN FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // DBARD EVEN BACKWARD
    cout << "Computing DBARD EVEN BACKWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double dbardEB;
    text = "listdbdcutEB[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            dbardEB = integration_dbardEB_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << dbardEB;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DBARD EVEN BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // DBARD ODD FORWARD
    cout << "Computing DBARD ODD FORWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double dbardOF;
    text = "listdbdcutOF[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            dbardOF = integration_dbardOF_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << dbardOF;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DBARD ODD FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    // DBARD ODD BACKWARD
    cout << "Computing DBARD ODD BACKWARD" << endl;
    begin = std::chrono::steady_clock::now();
    double dbardOB;
    text = "listdbdcutOB[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
    
    for (size_t PDF_set = 0; PDF_set <= 0; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            dbardOB = integration_dbardOB_y(PDF_set, Minv_table[i]);
                        
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << dbardOB;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
//     text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DBARD ODD BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    
    if (PDF_error_flag) {
        // UUBAR EVEN
        cout << "Computing UUBAR EVEN PDF ERROR" << endl;
        begin = std::chrono::steady_clock::now();
        
        double DeltaPDFuuE;
        
        text = "DeltaPDFuuE[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
        
        for (int i = 0; i < points_Minv; i++) {
            
            DeltaPDFuuE = 0;
            
            for (size_t PDF_set = 1; PDF_set <= nmem; PDF_set = PDF_set + 2) {
                
                DeltaPDFuuE += pow((integration_uubarEF_y(PDF_set, Minv_table[i]) + integration_uubarEB_y(PDF_set, Minv_table[i]) + integration_ubaruEF_y(PDF_set, Minv_table[i]) + integration_ubaruEB_y(PDF_set, Minv_table[i]))-(integration_uubarEF_y(PDF_set+1, Minv_table[i]) + integration_uubarEB_y(PDF_set+1, Minv_table[i]) + integration_ubaruEF_y(PDF_set+1, Minv_table[i]) + integration_ubaruEB_y(PDF_set+1, Minv_table[i])),2);
                
            }
                    
            DeltaPDFuuE = (0.5)*sqrt(DeltaPDFuuE);
            
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << DeltaPDFuuE;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
            }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
        
        text = text.substr(0, text.size()-3);
    //     text = text + "}";
        text = text + "\n\n";
        out << text;
        
        end = std::chrono::steady_clock::now();
        cout << "UUBAR EVEN PDF ERROR completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
        
        
//         // UUBAR ODD
//         cout << "Computing UUBAR ODD PDF ERROR" << endl;
//         begin = std::chrono::steady_clock::now();
//         
//         double DeltaPDFuuO;
//         
//         text = "DeltaPDFuuO[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             DeltaPDFuuO = 0;
//             
//             for (size_t PDF_set = 1; PDF_set <= nmem; PDF_set = PDF_set + 2) {
//                 
//                 DeltaPDFuuO += pow((integration_uubarOF_y(PDF_set, Minv_table[i]) + integration_uubarOB_y(PDF_set, Minv_table[i]) + integration_ubaruOF_y(PDF_set, Minv_table[i]) + integration_ubaruOB_y(PDF_set, Minv_table[i]))-(integration_uubarOF_y(PDF_set+1, Minv_table[i]) + integration_uubarOB_y(PDF_set+1, Minv_table[i]) + integration_ubaruOF_y(PDF_set+1, Minv_table[i]) + integration_ubaruOB_y(PDF_set+1, Minv_table[i])),2);
//                 
//             }
//                     
//             DeltaPDFuuO = (0.5)*sqrt(DeltaPDFuuO);
//             
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << DeltaPDFuuO;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//             }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//         
//         text = text.substr(0, text.size()-3);
//     //     text = text + "}";
//         text = text + "\n\n";
//         out << text;
//         
//         end = std::chrono::steady_clock::now();
//         cout << "UUBAR ODD PDF ERROR completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
        
        
        // DDBAR EVEN
        cout << "Computing DDBAR EVEN PDF ERROR" << endl;
        begin = std::chrono::steady_clock::now();
        
        double DeltaPDFddE;
        
        text = "DeltaPDFddE[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
        
        for (int i = 0; i < points_Minv; i++) {
            
            DeltaPDFddE = 0;
            
            for (size_t PDF_set = 1; PDF_set <= nmem; PDF_set = PDF_set + 2) {
                
                DeltaPDFddE += pow((integration_ddbarEF_y(PDF_set, Minv_table[i]) + integration_ddbarEB_y(PDF_set, Minv_table[i]) + integration_dbardEF_y(PDF_set, Minv_table[i]) + integration_dbardEB_y(PDF_set, Minv_table[i]))-(integration_ddbarEF_y(PDF_set+1, Minv_table[i]) + integration_ddbarEB_y(PDF_set+1, Minv_table[i]) + integration_dbardEF_y(PDF_set+1, Minv_table[i]) + integration_dbardEB_y(PDF_set+1, Minv_table[i])),2);
                
            }
                    
            DeltaPDFddE = (0.5)*sqrt(DeltaPDFddE);
            
            stringstream Minv_string;
            Minv_string << fixed << setprecision(3) << Minv_table[i];
            string Minv_s = Minv_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << DeltaPDFddE;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + Minv_s + ", " + data_s + "}, ";
            }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
        
        text = text.substr(0, text.size()-3);
    //     text = text + "}";
        text = text + "\n\n";
        out << text;
        
        end = std::chrono::steady_clock::now();
        cout << "DDBAR EVEN PDF ERROR completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
        
        
//         // DDBAR ODD
//         cout << "Computing DDBAR ODD PDF ERROR" << endl;
//         begin = std::chrono::steady_clock::now();
//         
//         double DeltaPDFddO;
//         
//         text = "DeltaPDFddO[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {";
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             DeltaPDFddO = 0;
//             
//             for (size_t PDF_set = 1; PDF_set <= nmem; PDF_set = PDF_set + 2) {
//                 
//                 DeltaPDFddO += pow((integration_ddbarOF_y(PDF_set, Minv_table[i]) + integration_ddbarOB_y(PDF_set, Minv_table[i]) + integration_dbardOF_y(PDF_set, Minv_table[i]) + integration_dbardOB_y(PDF_set, Minv_table[i]))-(integration_ddbarOF_y(PDF_set+1, Minv_table[i]) + integration_ddbarOB_y(PDF_set+1, Minv_table[i]) + integration_dbardOF_y(PDF_set+1, Minv_table[i]) + integration_dbardOB_y(PDF_set+1, Minv_table[i])),2);
//                 
//             }
//                     
//             DeltaPDFddO = (0.5)*sqrt(DeltaPDFddO);
//             
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << DeltaPDFddO;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//             }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//         
//         text = text.substr(0, text.size()-3);
//     //     text = text + "}";
//         text = text + "\n\n";
//         out << text;
//         
//         end = std::chrono::steady_clock::now();
//         cout << "DDBAR ODD PDF ERROR completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    }
    
    
    
    // Clear memory
    gsl_integration_workspace_free (w);

    out.close();
    
    return 0;
}
