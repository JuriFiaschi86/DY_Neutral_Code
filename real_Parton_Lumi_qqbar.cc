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

// set collider energy and luminosity
#define energy 13000

// PDF set and grid
#define setname "CT18NNLO68cl"
// #define setname "CT18NNLO68cl_AFBonAW_300"
// #define setname "CT18NNLO68cl_AFBonAW_3000"

const LHAPDF::PDFSet set(setname);
const size_t nmem = set.size()-1;
// const size_t nmem = 0;

const vector<LHAPDF::PDF*> pdf = set.mkPDFs();

// Invariant mass range and number of points
const double MX_min = 100;
const double MX_max = energy;
const int points_MX = 1200;
const double step = (MX_max - MX_min) / points_MX;

// Integration parameters
const int alloc_space = points_MX;
gsl_integration_workspace * w = gsl_integration_workspace_alloc(alloc_space);
// Integration option
const int key = 6;
// Error of the integration
const double epsabs = 0;
const double epsrel = 1e-2;

struct param_struct {
    double MX;
    int PDF_set;
};

////U-UBAR parton luminosity
double uu_funct (double x, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double MX = (p->MX);
    int PDF_set = (p->PDF_set);
    
    // Partonic cross section parameters
    double Q = MX;
    double tau = pow(MX/energy,2);
    
    // Partons PDFs uubar
    double f1u = (pdf[PDF_set]->xfxQ(2, x, Q))/x;
    double f2ubar = (pdf[PDF_set]->xfxQ(-2, tau/x, Q))/(tau/x);
    double f1c = (pdf[PDF_set]->xfxQ(4, x, Q))/x;
    double f2cbar = (pdf[PDF_set]->xfxQ(-4, tau/x, Q))/(tau/x);
    
    // Partons PDFs ubaru
    double f1ubar = (pdf[PDF_set]->xfxQ(-2, x, Q))/x;
    double f2u = (pdf[PDF_set]->xfxQ(2, tau/x, Q))/(tau/x);
    double f1cbar = (pdf[PDF_set]->xfxQ(-4, x, Q))/x;    
    double f2c = (pdf[PDF_set]->xfxQ(4, tau/x, Q))/(tau/x);
    
    // PDF combinations
    double uubar_PDF = f1u*f2ubar + f1c*f2cbar;
    double ubaru_PDF = f1ubar*f2u + f1cbar*f2c;
    
    return (uubar_PDF + ubaru_PDF)/x;
}

////U-UBAR parton luminosity integration
double integration_uu (int PDF_set, double MX) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &uu_funct;
    F.params = &p;
    
    p.MX = MX;
    p.PDF_set = PDF_set;
    
    double tau = pow(MX/energy,2);
    
    double inf = tau;
    double sup = 1.0;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error);
    
    return (2*MX/pow(energy,2))*result;
}

////D-DBAR parton luminosity
double dd_funct (double x, void * params) {
    
    struct param_struct * p = (struct param_struct *)params;
    
    // Pass the invariant mass as parameter
    double MX = (p->MX);
    int PDF_set = (p->PDF_set);
    
    // Partonic cross section parameters
    double Q = MX;
    double tau = pow(MX/energy,2);
    
    
    // Partons PDFs ddbar
    double f1d = (pdf[PDF_set]->xfxQ(1, x, Q))/x;
    double f2dbar = (pdf[PDF_set]->xfxQ(-1, tau/x, Q))/(tau/x);
    double f1s = (pdf[PDF_set]->xfxQ(3, x, Q))/x;
    double f2sbar = (pdf[PDF_set]->xfxQ(-3, tau/x, Q))/(tau/x);
    double f1b = (pdf[PDF_set]->xfxQ(5, x, Q))/x;
    double f2bbar = (pdf[PDF_set]->xfxQ(-5, tau/x, Q))/(tau/x);
    
    // Partons PDFs dbard
    double f1dbar = (pdf[PDF_set]->xfxQ(-1, x, Q))/x;
    double f2d = (pdf[PDF_set]->xfxQ(1, tau/x, Q))/(tau/x);
    double f1sbar = (pdf[PDF_set]->xfxQ(-3, x, Q))/x;
    double f2s = (pdf[PDF_set]->xfxQ(3, tau/x, Q))/(tau/x);
    double f1bbar = (pdf[PDF_set]->xfxQ(-5, x, Q))/x;
    double f2b = (pdf[PDF_set]->xfxQ(5, tau/x, Q))/(tau/x);
    
    // PDF combinations    
    double ddbar_PDF = f1d*f2dbar + f1s*f2sbar + f1b*f2bbar;
    double dbard_PDF = f1dbar*f2d + f1sbar*f2s + f1bbar*f2b;
    
    return (ddbar_PDF + dbard_PDF)/x;
}

////D-DBAR parton luminosity integration
double integration_dd (int PDF_set, double MX) {
    
    double result, error;
    struct param_struct p;
    
    gsl_function F;
    F.function = &dd_funct;
    F.params = &p;
    
    p.MX = MX;
    p.PDF_set = PDF_set;
    
    double tau = pow(MX/energy,2);
    
    double inf = tau;
    double sup = 1.0;
    
    gsl_integration_qag (&F, inf, sup, epsabs, epsrel, alloc_space, key, w, &result, &error);
    
    return (2*MX/pow(energy,2))*result;
}

int main(int argc, char** argv)
{
 
    // Time counter
    std::chrono::steady_clock::time_point begin, end;
    
    // Output file name
    char filename[80];
    strcpy(filename, "Parton_lumi_");
    strcat(filename, setname);
    strcat(filename, ".dat");
    
    ofstream out(filename);    
    
    vector<double> MX_table;
    // Invariant mass table
    MX_table.push_back(MX_min);
    for (int i = 1; i < points_MX; i++) {
        MX_table.push_back(MX_table[i-1]+step);
    }
    
    stringstream energy_string;
    energy_string << fixed << setprecision(0) << energy;
    string energy_s = energy_string.str();
    
    string text;
    
    // UU parton luminosity
    cout << "Computing UU parton luminosity" << endl;
    begin = std::chrono::steady_clock::now();
    double uu;
    text = "listuu[" + energy_s + "] = {{";
    
    for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
        
        for (int i = 0; i < points_MX; i++) {
            
            uu = integration_uu(PDF_set, MX_table[i]);
                        
            stringstream MX_string;
            MX_string << fixed << setprecision(3) << MX_table[i];
            string MX_s = MX_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << uu;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + MX_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
    text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "UU parton luminosity completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    // DD parton luminosity
    cout << "Computing DD parton luminosity" << endl;
    begin = std::chrono::steady_clock::now();
    double dd;
    text = "listdd[" + energy_s + "] = {{";
    
    for (size_t PDF_set = 0; PDF_set <= nmem; PDF_set++) {
        
        for (int i = 0; i < points_MX; i++) {
            
            dd = integration_dd(PDF_set, MX_table[i]);
                        
            stringstream MX_string;
            MX_string << fixed << setprecision(3) << MX_table[i];
            string MX_s = MX_string.str();
            
            stringstream data_string;
            data_string << fixed << scientific << setprecision(16) << dd;
            string data_s = data_string.str();
            data_s = regex_replace(data_s, regex("e"), "*^");
            
            text = text + "{" + MX_s + ", " + data_s + "}, ";
        }
        text = text.substr(0, text.size()-2);
        text = text + "}, {";
    }
    
    text = text.substr(0, text.size()-3);
    text = text + "}";
    text = text + "\n\n";
    out << text;
    
    end = std::chrono::steady_clock::now();
    cout << "DD parton luminosity completed in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " sec." << endl;
    
    // Clear memory
    gsl_integration_workspace_free (w);

    out.close();
    
    return 0;
}
