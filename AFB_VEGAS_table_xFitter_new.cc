#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <time.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "LHAPDF/LHAPDF.h"

using namespace LHAPDF;
using namespace std;

#define verbosity 0

// // SM input (v0 Mathematica)
// // Constants
// #define PI 3.14159265
// #define GeVtofb 0.38937966e+12
// 
// // SM parameters
// #define MZ 91.19
// #define GammaZ 2.5
// #define alphaEM 7.8125e-3
// #define stheta2W 0.232035

// SM input (v1 xFitter)
// Constants
#define PI 3.14159265
#define GeVtofb 0.389379338e+12

// SM parameters
#define MZ 91.1876
#define GammaZ 2.4952
#define alphaEM 7.29735e-3
#define stheta2W 0.23127

// //// SM input (v2)
// // Constants
// #define PI 3.14159265
// #define GeVtofb 0.38937966e+12
// 
// // SM parameters
// #define MZ 91.19
// #define GammaZ 2.5
// #define alphaEM 7.8125e-3
// #define stheta2W 0.232035

// // SM input (v3 MCFM)
// #define PI 3.14159265
// #define GeVtofb 0.38937966e+12
// 
// #define MZ 91.1876
// #define GammaZ 2.4952
// #define alphaEM 7.756247e-3
// #define stheta2W 0.2228972


// set collider energy and luminosity
#define energy 13000
#define lum 3000

// define acceptance cuts
#define eta_cut 2.5
#define pT_cut 20

// define rapidity cuts
#define y_min 0.0
#define y_max 0.0 // y_max = 0 means no upper cut

// PDF set and grid
// #define setname "CT18NNLO"
// #define setname "NNPDF30_nlo_as_0118_hessian"
#define setname "NNPDF40_nnlo_as_01180"
#define PDF_set 0

// PDFs section
// Initialize PDFs
void initLHAPDF();
const vector<LHAPDF::PDF*> pdf = mkPDFs(setname);


// Setting of the integration
const int dim_integration = 2; // Integration in yreduced and Minv
// Integration extremes
double yreducedmin = 0;
double yreducedmax = 1;

// Integration number of calls
size_t calls = 10000;
#define max_iter 10


double *propagators (double Minv)
{
    const double e = sqrt(4*PI*alphaEM);
    const double gsm = e/(sqrt(stheta2W));

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
//     couplings_Z[0] = (1.0/2.0)*gsm*(1.0/6.0)*(3*cos(smangle)+8*sin(smangle));
//     couplings_Z[1] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
//     couplings_Z[2] = (1.0/2.0)*gsm*(1.0/6.0)*(-3*cos(smangle)-4*sin(smangle));
//     couplings_Z[3] = (1.0/2.0)*gsm*(-cos(smangle)/2.0);
//     couplings_Z[4] = (1.0/2.0)*gsm*((-cos(smangle)/2.0)+(-2*sin(smangle)));
//     couplings_Z[5] = (1.0/2.0)*gsm*(-cos(smangle)/2.0);
//     couplings_Z[6] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
//     couplings_Z[7] = (1.0/2.0)*gsm*(cos(smangle)/2.0);
    
    
    couplings_Z[0] = (1.0/2.0)*gsm/sqrt(1 - stheta2W)*(1.0/2.0 - 2*(2.0/3.0)*stheta2W);
    couplings_Z[1] = (1.0/2.0)*gsm/sqrt(1 - stheta2W)*(1.0/2.0);
    couplings_Z[2] = (1.0/2.0)*gsm/sqrt(1 - stheta2W)*(-1.0/2.0 - 2*(-1.0/3.0)*stheta2W);
    couplings_Z[3] = (1.0/2.0)*gsm/sqrt(1 - stheta2W)*(-1.0/2.0);
    couplings_Z[4] = (1.0/2.0)*gsm/sqrt(1 - stheta2W)*(-1.0/2.0 - 2*(-1.0)*stheta2W);
    couplings_Z[5] = (1.0/2.0)*gsm/sqrt(1 - stheta2W)*(-1.0/2.0);
    couplings_Z[6] = (1.0/2.0)*gsm/sqrt(1 - stheta2W)*(1.0/2.0 - 2*(0.0)*stheta2W);
    couplings_Z[7] = (1.0/2.0)*gsm/sqrt(1 - stheta2W)*(1.0/2.0);

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
double uubarEF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));    
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////UUBAR EVEN FORWARD Integration
double integration_uubarEF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};
    double uubarEF, error_uubarEF;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_uubarEF = { &uubarEF_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_uubarEF, integration_inf, integration_sup, dim_integration, calls, r, s, &uubarEF, &error_uubarEF);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_uubarEF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &uubarEF, &error_uubarEF);
        ii++;
        if(ii = max_iter) break;
    }
    while ((fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5));
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*uubarEF;
}

////UUBAR EVEN BACKWARD Matrix element
double uubarEB_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////UUBAR EVEN BACKWARD Integration
double integration_uubarEB (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double uubarEB, error_uubarEB;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_uubarEB = { &uubarEB_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_uubarEB, integration_inf, integration_sup, dim_integration, calls, r, s, &uubarEB, &error_uubarEB);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_uubarEB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &uubarEB, &error_uubarEB);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*uubarEB;
}

////UUBAR ODD FORWARD Matrix element
double uubarOF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////UUBAR ODD FORWARD Integration
double integration_uubarOF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double uubarOF, error_uubarOF;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_uubarOF = { &uubarOF_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_uubarOF, integration_inf, integration_sup, dim_integration, calls, r, s, &uubarOF, &error_uubarOF);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_uubarOF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &uubarOF, &error_uubarOF);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*uubarOF;
}

////UUBAR ODD BACKWARD Matrix element
double uubarOB_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////UUBAR ODD BACKWARD Integration
double integration_uubarOB (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double uubarOB, error_uubarOB;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_uubarOB = { &uubarOB_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_uubarOB, integration_inf, integration_sup, dim_integration, calls, r, s, &uubarOB, &error_uubarOB);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_uubarOB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &uubarOB, &error_uubarOB);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*uubarOB;
}

////UBARU EVEN FORWARD Matrix element
double ubaruEF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////UBARU EVEN FORWARD Integration
double integration_ubaruEF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double ubaruEF, error_ubaruEF;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_ubaruEF = { &ubaruEF_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_ubaruEF, integration_inf, integration_sup, dim_integration, calls, r, s, &ubaruEF, &error_ubaruEF);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_ubaruEF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ubaruEF, &error_ubaruEF);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*ubaruEF;
}

////UBARU EVEN BACKWARD Matrix element
double ubaruEB_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////UBARU EVEN BACKWARD Integration
double integration_ubaruEB (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double ubaruEB, error_ubaruEB;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_ubaruEB = { &ubaruEB_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_ubaruEB, integration_inf, integration_sup, dim_integration, calls, r, s, &ubaruEB, &error_ubaruEB);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_ubaruEB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ubaruEB, &error_ubaruEB);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*ubaruEB;
}

////UBARU ODD FORWARD Matrix element
double ubaruOF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////UBARU ODD FORWARD Integration
double integration_ubaruOF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double ubaruOF, error_ubaruOF;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_ubaruOF = { &ubaruOF_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_ubaruOF, integration_inf, integration_sup, dim_integration, calls, r, s, &ubaruOF, &error_ubaruOF);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_ubaruOF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ubaruOF, &error_ubaruOF);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*ubaruOF;
}

////UBARU ODD BACKWARD Matrix element
double ubaruOB_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////UBARU ODD BACKWARD Integration
double integration_ubaruOB (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double ubaruOB, error_ubaruOB;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_ubaruOB = { &ubaruOB_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_ubaruOB, integration_inf, integration_sup, dim_integration, calls, r, s, &ubaruOB, &error_ubaruOB);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_ubaruOB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ubaruOB, &error_ubaruOB);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*ubaruOB;
}

////DDBAR EVEN FORWARD Matrix element
double ddbarEF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////DDBAR EVEN FORWARD Integration
double integration_ddbarEF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double ddbarEF, error_ddbarEF;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_ddbarEF = { &ddbarEF_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_ddbarEF, integration_inf, integration_sup, dim_integration, calls, r, s, &ddbarEF, &error_ddbarEF);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_ddbarEF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ddbarEF, &error_ddbarEF);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*ddbarEF;
}

////DDBAR EVEN BACKWARD Matrix element
double ddbarEB_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////DDBAR EVEN BACKWARD Integration
double integration_ddbarEB (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double ddbarEB, error_ddbarEB;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_ddbarEB = { &ddbarEB_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_ddbarEB, integration_inf, integration_sup, dim_integration, calls, r, s, &ddbarEB, &error_ddbarEB);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_ddbarEB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ddbarEB, &error_ddbarEB);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*ddbarEB;
}

////DDBAR ODD FORWARD Matrix element
double ddbarOF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////DDBAR ODD FORWARD Integration
double integration_ddbarOF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double ddbarOF, error_ddbarOF;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_ddbarOF = { &ddbarOF_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_ddbarOF, integration_inf, integration_sup, dim_integration, calls, r, s, &ddbarOF, &error_ddbarOF);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_ddbarOF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ddbarOF, &error_ddbarOF);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*ddbarOF;
}

////DDBAR ODD BACKWARD Matrix element
double ddbarOB_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////DDBAR ODD BACKWARD Integration
double integration_ddbarOB (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double ddbarOB, error_ddbarOB;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_ddbarOB = { &ddbarOB_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_ddbarOB, integration_inf, integration_sup, dim_integration, calls, r, s, &ddbarOB, &error_ddbarOB);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_ddbarOB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &ddbarOB, &error_ddbarOB);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*ddbarOB;
}


////DBARD EVEN FORWARD Matrix element
double dbardEF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////DBARD EVEN FORWARD Integration
double integration_dbardEF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double dbardEF, error_dbardEF;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_dbardEF = { &dbardEF_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_dbardEF, integration_inf, integration_sup, dim_integration, calls, r, s, &dbardEF, &error_dbardEF);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_dbardEF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &dbardEF, &error_dbardEF);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*dbardEF;
}

////DBARD EVEN BACKWARD Matrix element
double dbardEB_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////DBARD EVEN BACKWARD Integration
double integration_dbardEB (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double dbardEB, error_dbardEB;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_dbardEB = { &dbardEB_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_dbardEB, integration_inf, integration_sup, dim_integration, calls, r, s, &dbardEB, &error_dbardEB);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_dbardEB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &dbardEB, &error_dbardEB);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*dbardEB;
}

////DBARD ODD FORWARD Matrix element
double dbardOF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////DBARD ODD FORWARD Integration
double integration_dbardOF (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double dbardOF, error_dbardOF;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_dbardOF = { &dbardOF_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_dbardOF, integration_inf, integration_sup, dim_integration, calls, r, s, &dbardOF, &error_dbardOF);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_dbardOF, integration_inf, integration_sup, dim_integration, calls/5, r, s, &dbardOF, &error_dbardOF);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*dbardOF;
}

////DBARD ODD BACKWARD Matrix element
double dbardOB_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];
    double Minv = entries[1];

    // Partonic cross section parameters
    double Q = Minv;
    double z = pow(Minv,2)/pow(energy,2);
    double y = -(1.0/2.0)*log(z)*(yreduced);
    double x1 = sqrt(z)*exp(y);
    double x2 = sqrt(z)*exp(-y);
    double dsigma_temp = pow(Minv,2)/(96*PI);
    double dsigma = GeVtofb*dsigma_temp*(2*Minv/pow(energy,2))*(-(1.0/2.0)*log(z));
    
    if (y < y_min or (y > y_max and y_max != 0)) {
        dsigma = 0;
    }

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

////DBARD ODD BACKWARD Integration
double integration_dbardOB (double Minv_inf, double Minv_sup)
{
    double integration_inf[2] = {yreducedmin, Minv_inf};
    double integration_sup[2] = {yreducedmax, Minv_sup};

    double dbardOB, error_dbardOB;
    // Initialization of the integration (quite a black box)
    gsl_monte_function Integrate_dbardOB = { &dbardOB_funct, dim_integration, 0 };
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    // Integration
    {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim_integration);
    gsl_monte_vegas_integrate (&Integrate_dbardOB, integration_inf, integration_sup, dim_integration, calls, r, s, &dbardOB, &error_dbardOB);
    int ii = 0;
    do
    {
        gsl_monte_vegas_integrate (&Integrate_dbardOB, integration_inf, integration_sup, dim_integration, calls/5, r, s, &dbardOB, &error_dbardOB);
        ii++;
        if(ii = max_iter) break;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free (s);
    }

    gsl_rng_free (r);

    return 2*dbardOB;
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


    double *results = new double[6];

    results[0] = XS;
    results[1] = Epsilon_XS;
    results[2] = AFB;
    results[3] = Delta_AFB;
    
    results[4] = Forward;
    results[5] = Backward;

    return results;
}

int main(int argc, char** argv)
{
    double Minv_min, Minv_max;
    double resolution;
    int size;
    
    string path_APPLgrids = "/home/juri/xfitter_before_ceres_merge/xfitter-master_before_PionCeres_merge/datafiles/AFB/AFB_13TeV/";
    
    Minv_min = 80.0;
    Minv_max = 200.0;
    
    resolution = 2.5;
    
    size = int((Minv_max - Minv_min)/resolution);
    
    // check on the rapidity cut
    if (y_min >= eta_cut) {
        printf("\nThe chosen lower rapidity cut is higher than pseudorapidity cut.\n\n");
        return 0;
    }
    if (y_min / log(energy/Minv_max) > 1) {
        printf("\nThe chosen lower rapidity cut is too high in this invariant mass range.\n\n");
        return 0;
    }
    
    double Minv_inf[size];
    double Minv_sup[size];
    
    double AFB_values[size];
    double Delta_AFB_values[size];
    
    for (int i = 0; i < size; i++) {
        Minv_inf[i] = Minv_min + i*resolution;
        Minv_sup[i] = Minv_min + (i+1)*resolution;
        double *results = observables (Minv_inf[i], Minv_sup[i]);
        AFB_values[i] = results[2];
        Delta_AFB_values[i] = results[3];
    }
    
    // contruct the file name of the datafile    
    char const *setname_string;
    if (not(strcmp(setname, "CT14nnlo"))) {
        setname_string = "CT14nnlo";
    }
    else if (not(strcmp(setname, "CT18NNLO"))) {
        setname_string = "CT18nnlo";
    }
    else if (not(strcmp(setname, "ABMP16_5_nnlo"))) {
        setname_string = "ABMP16nnlo";
    }
    else if (not(strcmp(setname, "HERAPDF20_NNLO_EIG"))) {
        setname_string = "HERAPDF20nnlo";
    }
    else if (not(strcmp(setname, "MMHT2014nnlo68cl"))) {
        setname_string = "MMHT2014nnlo";
    }
    else if (not(strcmp(setname, "NNPDF31_nnlo_as_0118_hessian"))) {
        setname_string = "NNPDF31nnlo";
    }
    
    char l[10];
    sprintf(l, "%d", lum);
    char filename[30];
    strcpy(filename, "AFB_");
    strcat(filename, setname_string);
    strcat(filename, "_");
    strcat(filename, l);
    strcat(filename, ".dat");
    
    // get current date
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char date[11];
    sprintf(date, "%d-%d-%d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
    
    ofstream out(filename);
    
    string text = " !*  File produced by Juri Fiaschi\n !*  Created on \0";
    out << text;
    out << date;
    out << "\n\n";
    out << "&Data\n  Name = \"DY AFB\"\n  IndexDataset = 102\n  Reaction = \"DY AFB pp->l+l-\"\n\n  TheoryType = 'expression'\n  TermName = 'PY0','PY1','PY2','PY3','PY4','MY0','MY1','MY2','MY3','MY4'\n  TermType = 'reaction','reaction','reaction','reaction','reaction','reaction','reaction','reaction','reaction','reaction'\n  TermSource = 'APPLgrid','APPLgrid','APPLgrid','APPLgrid','APPLgrid','APPLgrid','APPLgrid','APPLgrid','APPLgrid','APPLgrid'\n";
    out << "  TermInfo =  'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_POS/aMCfast_obs_0.root," << path_APPLgrids << "CosTheta_POS_SecondBin/aMCfast_obs_0.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_POS/aMCfast_obs_1.root," << path_APPLgrids << "CosTheta_POS_SecondBin/aMCfast_obs_1.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_POS/aMCfast_obs_2.root," << path_APPLgrids << "CosTheta_POS_SecondBin/aMCfast_obs_2.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_POS/aMCfast_obs_3.root," << path_APPLgrids << "CosTheta_POS_SecondBin/aMCfast_obs_3.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_POS/aMCfast_obs_4.root," << path_APPLgrids << "CosTheta_POS_SecondBin/aMCfast_obs_4.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_NEG/aMCfast_obs_0.root," << path_APPLgrids << "CosTheta_NEG_SecondBin/aMCfast_obs_0.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_NEG/aMCfast_obs_1.root," << path_APPLgrids << "CosTheta_NEG_SecondBin/aMCfast_obs_1.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_NEG/aMCfast_obs_2.root," << path_APPLgrids << "CosTheta_NEG_SecondBin/aMCfast_obs_2.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_NEG/aMCfast_obs_3.root," << path_APPLgrids << "CosTheta_NEG_SecondBin/aMCfast_obs_3.root',\n";
    out << "              'Order=NLO:GridName=" << path_APPLgrids << "CosTheta_NEG/aMCfast_obs_4.root," << path_APPLgrids << "CosTheta_NEG_SecondBin/aMCfast_obs_4.root',\n";
    out << "  TheorExpr = '(PY0+PY1+PY2+PY3+PY4-MY0-MY1-MY2-MY3-MY4)/(PY0+PY1+PY2+PY3+PY4+MY0+MY1+MY2+MY3+MY4)'\n";
    out << "\n";
    out << "  NData =    " << size << "\n";
    out << "  NColumn =   4\n";
    out << "  ColumnType = 2*\"Bin\",\"Sigma\", 1*\"Error\"\n";
    out << "  ColumnName = \"Minv_min\", \"Minv_max\", \"Sigma\", \"stat\"\n";
    out << "&End\n";
    out << "&PlotDesc\n";
    out << "         PlotN = 1\n";
    out << "         PlotDefColumn = 'Minv_max'\n";
    out << "         PlotDefValue = " << Minv_min << ", " << Minv_max << "\n";
    out << "         PlotOptions(1)  = 'Experiment:pseudodata@ExtraLabel:pp #rightarrow Z; #sqrt{s} = 13 TeV@Title: @XTitle: M(ll) @YTitle: AFB '\n";
    out << "&End\n";
    out << "*        Minv_min           Minv_max           Sigma    stat\n";
    cout << "Minv_min\tMinv_max\tAFB*\t\t\u0394AFB*\n";
    for (int i = 0; i < size; i++) {
        out << fixed << setprecision(1) << Minv_inf[i] << "\t" << Minv_sup[i] << "\t" << scientific << setprecision(4) << AFB_values[i] << "\t" << Delta_AFB_values[i] << "\n";
        cout << fixed << setprecision(1) << Minv_inf[i] << "\t\t" << Minv_sup[i] << "\t\t" << scientific << setprecision(4) << AFB_values[i] << "\t" << Delta_AFB_values[i] << "\n";
    }
    
    out.close();
    
    cout << "\nCreated " << filename << " datafile for xFitter\n\n" << endl;
    
    return 0;
}
