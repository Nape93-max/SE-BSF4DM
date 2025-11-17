#include"../../Packages/SE_BSF/SE_BSF_header.h"
#include"../../Packages/SE_BSF/SE_BSF_functions.cpp"

/* In this file, the user has to supply some details of the model for BSF. The following information needs to be provided:
 1) The BSF scenario/limit
 2) The number of excited states to be considered
*/

int bsf_scenario = 0; //Flag for BSF (0 = no BSF, 1 = no transition limit, 2 = efficient transition limit, 3 = ionization equilibrium, 4 = Full matrix solution)
int num_excited_states = 0; //Number of n states included in the calculation (n = 0 -> no bound states, n = 1 -> ground state etc.)

// ************************** END USER DEFINITION *********************************************

int num_of_mediators; 
int* pdg_nums_mediators = getDarkSectorPDGs(&num_of_mediators); //Retrieve the number of dark sector particles and their PDG numbers

double BoundStateCoannihilation(double T){ /* This function can be modified by the user according to the theory at hand. 
It calculates the (co-)annihilation terms for bound state formation in the spirit of eqns. (2.3) - (2.10) of arXiv:2203.04326v2.
If left unchanged, it calculates the BSF in all possible coannihilation pairings.    */

    double bsf_xs = 0.; //Initialise final output 
    double xJH = Mcdm/T;

    int num_of_coann_terms = (num_of_mediators*(num_of_mediators + 1))/2; // N*(N+1)/2 different terms for N mediators
    
    double mediator_masses[num_of_mediators]; //Assign masses to the mediators
    double deltaJHarray[num_of_mediators]; //Values for delta for each mediator according to eq. (2.9) of arXiv:2203.04326v2
    double YeqContribs[num_of_mediators]; //Assign values for Yeq for mediator i according to eq. (2.8) of arXiv:2203.04326v2
    double YeqTilde = 0.; //Quantity defined in eq. (2.3) of arXiv:2203.04326v2

    for(int i=0; i<num_of_mediators; i++){ //Initializations and definitions
        char* name_temp = new char[strlen(pdg2name(pdg_nums_mediators[i])) + 1]; // +1 for null terminator
        if (name_temp != nullptr) {
            strcpy(name_temp, pdg2name(pdg_nums_mediators[i]));
        } 
        mediator_masses[i] = pMass(name_temp);
        deltaJHarray[i] = mediator_masses[i]/Mcdm - 1.;
        YeqContribs[i] = 0.;

        if(deltaJHarray[i] < 0.5){ //Take only sufficiently light d.o.f.s into account
            YeqContribs[i] = dofs_particle(name_temp)*pow(1. + deltaJHarray[i], 3./2.)*exp(-xJH*deltaJHarray[i]); //see eqns. (2.7 - 2.8) of arXiv:2203.04326v2
            YeqTilde += YeqContribs[i];
        }
        delete[] name_temp;
    }
    double YeqTildeSquared = YeqTilde*YeqTilde;

    for(int i=0; i<num_of_mediators; i++){ //This performs the coannihilation sum
        if(deltaJHarray[i] < 0.5){
            double bsf_part = BSF_XS_TA_ex(pdg_nums_mediators[i], pdg_nums_mediators[i], T, bsf_scenario, num_excited_states)*YeqContribs[i]*YeqContribs[i];
            bsf_xs += bsf_part; //This part counts only once!
        }
        for(int j=i+1; j<num_of_mediators; j++){
            if(deltaJHarray[i] < 0.5 && deltaJHarray[j] < 0.5){
                double bsf_part = BSF_XS_TA_ex(pdg_nums_mediators[i], pdg_nums_mediators[j], T, bsf_scenario, num_excited_states)*YeqContribs[i]*YeqContribs[j];
                bsf_xs += 2.*bsf_part; //This part counts twice!
            }
        }
    }

    return bsf_xs/YeqTildeSquared; 
}

// FREE the memory when done!
//free(pdg_nums_test);