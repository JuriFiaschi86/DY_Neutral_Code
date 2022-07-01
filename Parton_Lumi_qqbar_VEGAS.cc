#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <regex>
#include <chrono>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
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

// define rapidity cuts
#define y_min 0.0
#define y_max 0.0 // y_max = 0 means no upper cut

// PDF set and grid
#define setname "CT18NNLO68cl"
// #define PDF_set 0



const LHAPDF::PDFSet set(setname);
const size_t nmem = set.size()-1;

const vector<LHAPDF::PDF*> pdf = set.mkPDFs();

// Invariant mass range and number of points
// const double Minv_min = 40.001;
const double Minv_min = 2*pT_cut + 0.001;
const double Minv_max = energy;
const int points_Minv = 2;
const double step = (Minv_max - Minv_min) / points_Minv;
// const vector<double> Minv_table;

// Setting of the integration
const int dim_integration = 1; // Integration in yreduced and Minv
// Integration extremes
double yreducedmin = 0;
double yreducedmax = 1;

// Integration number of calls
size_t calls = 10000;
#define max_iter 10

struct param_struct {
    double Minv;
    int PDF_set;
};


////UUBAR EVEN FORWARD Matrix element
double uubarEF_funct (double *entries, size_t dim, void *params)
{
    (void)(dim); /* avoid unused parameter warnings */
    double yreduced = entries[0];

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

    return uubarEF;
}

////UUBAR EVEN FORWARD Integration
double integration_uubarEF (int PDF_set, double Minv)
{
    double integration_inf[1] = {yreducedmin};
    double integration_sup[1] = {yreducedmax};
    double uubarEF, error_uubarEF;
    
    struct param_struct p;
    p.Minv = Minv;
    p.PDF_set = PDF_set;
    
    // Initialization of the integration (quite a black box)
    
    gsl_monte_function Integrate_uubarEF = { &uubarEF_funct, dim_integration, &p };
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



int main(int argc, char** argv)
{
 
    // Time counter
    std::chrono::steady_clock::time_point begin, end;
    
    // Output file name
    char filename[40];
    strcpy(filename, "DY_");
    strcat(filename, setname);
    strcat(filename, "_eigen_VEGAS.dat");
    
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
    text = text + "listuubcutEFeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
    
    for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
        
        for (int i = 0; i < points_Minv; i++) {
            
            uubarEF = integration_uubarEF(PDF_set, Minv_table[i]);
                        
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
    text = text + "}";
    text = text + "\n\n";
    
    end = std::chrono::steady_clock::now();
    cout << "UUBAR EVEN FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;

    
//     // UUBAR EVEN BACKWARD
//     cout << "Computing UUBAR EVEN BACKWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double uubarEB;
//     text = text + "listuubcutEBeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             uubarEB = integration_uubarEB_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << uubarEB;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "UUBAR EVEN BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // UUBAR ODD FORWARD
//     cout << "Computing UUBAR ODD FORWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double uubarOF;
//     text = text + "listuubcutOFeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             uubarOF = integration_uubarOF_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << uubarOF;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "UUBAR ODD FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // UUBAR ODD BACKWARD
//     cout << "Computing UUBAR ODD BACKWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double uubarOB;
//     text = text + "listuubcutOBeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             uubarOB = integration_uubarOB_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << uubarOB;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "UUBAR ODD BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // UBARU EVEN FORWARD
//     cout << "Computing UBARU EVEN FORWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double ubaruEF;
//     text = text + "listubucutEFeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             ubaruEF = integration_ubaruEF_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << ubaruEF;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "UBARU EVEN FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // UBARU EVEN BACKWARD
//     cout << "Computing UBARU EVEN BACKWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double ubaruEB;
//     text = text + "listubucutEBeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             ubaruEB = integration_ubaruEB_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << ubaruEB;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "UBARU EVEN BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // UBARU ODD FORWARD
//     cout << "Computing UBARU ODD FORWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double ubaruOF;
//     text = text + "listubucutOFeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             ubaruOF = integration_ubaruOF_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << ubaruOF;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "UBARU ODD FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // UBARU ODD BACKWARD
//     cout << "Computing UBARU ODD BACKWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double ubaruOB;
//     text = text + "listubucutOBeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             ubaruOB = integration_ubaruOB_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << ubaruOB;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "UBARU ODD BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // DDBAR EVEN FORWARD
//     cout << "Computing DDBAR EVEN FORWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double ddbarEF;
//     text = text + "listddbcutEFeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             ddbarEF = integration_ddbarEF_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << ddbarEF;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "DDBAR EVEN FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // DDBAR EVEN BACKWARD
//     cout << "Computing DDBAR EVEN BACKWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double ddbarEB;
//     text = text + "listddbcutEBeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             ddbarEB = integration_ddbarEB_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << ddbarEB;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "DDBAR EVEN BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // DDBAR ODD FORWARD
//     cout << "Computing DDBAR ODD FORWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double ddbarOF;
//     text = text + "listddbcutOFeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             ddbarOF = integration_ddbarOF_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << ddbarOF;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "DDBAR ODD FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // DDBAR ODD BACKWARD
//     cout << "Computing DDBAR ODD BACKWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double ddbarOB;
//     text = text + "listddbcutOBeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             ddbarOB = integration_ddbarOB_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << ddbarOB;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "DDBAR ODD BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // DBARD EVEN FORWARD
//     cout << "Computing DBARD EVEN FORWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double dbardEF;
//     text = text + "listdbdcutEFeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             dbardEF = integration_dbardEF_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << dbardEF;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "DBARD EVEN FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // DBARD EVEN BACKWARD
//     cout << "Computing DBARD EVEN BACKWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double dbardEB;
//     text = text + "listdbdcutEBeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             dbardEB = integration_dbardEB_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << dbardEB;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "DBARD EVEN BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // DBARD ODD FORWARD
//     cout << "Computing DBARD ODD FORWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double dbardOF;
//     text = text + "listdbdcutOFeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             dbardOF = integration_dbardOF_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << dbardOF;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "DBARD ODD FORWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
//     
//     
//     // DBARD ODD BACKWARD
//     cout << "Computing DBARD ODD BACKWARD" << endl;
//     begin = std::chrono::steady_clock::now();
//     double dbardOB;
//     text = text + "listdbdcutOBeigen[" + eta_cut_s + ", " + pT_cut_s + "][" + energy_s + "] = {{";
//     
//     for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
//         
//         for (int i = 0; i < points_Minv; i++) {
//             
//             dbardOB = integration_dbardOB_y(PDF_set, Minv_table[i]);
//                         
//             stringstream Minv_string;
//             Minv_string << fixed << setprecision(3) << Minv_table[i];
//             string Minv_s = Minv_string.str();
//             
//             stringstream data_string;
//             data_string << fixed << scientific << setprecision(16) << dbardOB;
//             string data_s = data_string.str();
//             data_s = regex_replace(data_s, regex("e"), "*^");
//             
//             text = text + "{" + Minv_s + ", " + data_s + "}, ";
//         }
//         text = text.substr(0, text.size()-2);
//         text = text + "}, {";
//     }
//     
//     text = text.substr(0, text.size()-3);
//     text = text + "}";
//     text = text + "\n\n";
//     
//     end = std::chrono::steady_clock::now();
//     cout << "DBARD ODD BACKWARD completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;

    // Clear memory
//     gsl_integration_workspace_free (w);

    // Save to text
    out << text;
    out.close();
    
    return 0;
}
