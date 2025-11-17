
#ifndef SE_BSF_H
#define SE_BSF_H

#include"../../include/micromegas.h"
#include"../../include/micromegas_aux.h"
#include<array>

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

/* These 4 packages are required by ChatGPT for the DD file. */
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

typedef struct
{
    long n1;
    long n2;
    int n;
    int l;
    double T;
    double as_ann;
    double ag_B;
    double m1;
    double m2;
    double mu;
} Parameters;

typedef struct 
{
    int n;
    int l;
    int np;
    int lp;
    double mu;
    double alpha;
    double alpha_p;
} TransitionParameters;

int* getDarkSectorPDGs(int *num);
double smandel(double v, double m1, double m2);
double plowv(double v, double m1, double m2);
double smandelp(double Pcm, double m1, double m2);
double vRel(double Pcm, double m1, double m2);
double SommerfeldCoulomb(double zeta);
double BSF_XS_A(double T);
double BSF_XS_S(double T);
double BoundStateCoannihilation(double T);
double Ebin(double mu, double alpha, int n_BSF);
double GammaDecayAnalytic(Parameters pars);
double BS_mass(double m1, double m2, double Ebin);
double BS_dofs(int sb, int cb, int l);
double casimir2(int color);
double alphaQCD_binding(double fac);
double alpha_group(int r1, int r2, int r, double alpha_s);
double Pnl2_prefac_squared(double zeta_s, double zeta_bt, int n, int l);
double IR_bracket_part_high(double zeta_s, double zeta_bt, std::array<double, 5> Rnl, int l);
double IR_bracket_part_low(double zeta_s, double zeta_bt, std::array<double, 5> Rnl, int l);
double Gamma_Ion_Milne_relation(double sigma_v_averaged, Parameters pars);
double BSF_XS_TA_ex(long n1, long n2, double T, int scenario, int n_BSF);
double dofs_particle(char* prtcl_name);
double GammaIonGenMilne(double sigma_v_sum, Parameters pars);
double Gamma_Decay_effective_ET(Parameters pars);
double GammaTrans(TransitionParameters pars);
double hypergeometric_prior(double a, double b, double c, double z);
double* MatrixInversion(int dim, long double *MInput);
double Sigma_v_Ionization_Equilibrium(Parameters pars);
std::array<double, 5> Rnl_recursion(int n, int l, double zeta_bt, double zeta_s);
double bsf_cs_integrand_excited_states(double vrel, Parameters pars);
std::array<int, 2>index_conversion(int i);

#endif
