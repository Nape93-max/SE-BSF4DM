#include"../../Packages/SE_BSF/SE_BSF_header.h"

int somm_flag = 0;    // Flag for Sommerfeld effect (0 = not active, 1 = active, else = user defined)

double SommerfeldFactor_BSMmodel(double alphaQCD, double vrel, int c1, int c2, int c3, int c4, int s1, int s2, int s3, int s4, long n1, long n2, long n3, long n4)
{ /* This function calculates the Sommerfeld factor for a BSM model
     according to the color decomposition implemented by the user.
     alphaQCD is the value of the strong coupling at the appropriate scale, vrel the relative velocity,
     ci and ni; i = 1,2,3,4; are the dimensions of the SU(N) representation and the names of the particles in the process, respectively.
     Finally, the parameter cp indicates the s- and p-wave parts. cp = 0 (1) returns the factor for the s-wave (p-wave) cross section.
     
     Basically, the user has to has to adapt the "if...else" block according to the correct color decomposition. */

  double cfac1=4./3.; double cfac8=-1./6.; double cfac3=2./3.; double cfac6=-1./3.; //cfac is the coefficient of the coupling in the ARGUMENT of the Sommerfeld factor
  double kQfac1=0.; double kQfac8=0.; double kQfac3=0.; double kQfac6=0.; //kQfac is the coefficient IN FRONT OF the Sommerfeld factor coming from the color decomposition
  double zeta1=0.; double zeta3=0.; double zeta8=0.; double zeta6=0.; // zeta = alpha_group/vrel

// **************** BEGIN IF BLOCK OF COLOR DECOMPOSITION (CAN BE MODIFIED OPTIONALLY) ****************************************

  if((c1==3&&c2==-3)||(c1==-3&&c2==3)){ // Y Ybar process
    
    if((c3!=8&&c4==8)||(c3==8&&c4!=8)){ //g + Z/\gamma. This is purely adjoint for all partial waves. 
      kQfac8=1.;
      } 
    
    if(c3==1&&c4==1){ // most frequent case: Both final states are colour singlets
      kQfac1=1.; //for all partial waves 
      }

    if(c3==8&&c4==8){ //gg final state
      kQfac1=2./7.; kQfac8=5./7.; //gg channel
      }
    
    if((c3==3&&c4==-3)||(c3==-3&&c4==3)){ //q qbar final state
      /* This is the tricky channel, as elaborated in our publication.
      Interference terms in the color decomposition are neglected   */

      if((s1==1 && s2==1)&&(n3 == -n4)){ 
        // fermionic mediators and identical quark flavors in the final state
        kQfac8=1.; // for all partial waves in the case \lambda << g_s
      }
      else{
        // different quark flavors in the final state or scalar mediators
        kQfac1=1./9.; kQfac8=8./9.; //for all partial waves in the case \lambda << g_s
        }
      }
    }

  if((c1==3&&c2==3)||(c1==-3&&c2==-3)){
    // q_i q_j or q_i_bar, q_j_bar final state

    if((s1==0 && s2==0)&&(n3==n4)) { 
      // Scalar mediators and equal quark flavors in the final state
      kQfac6=1.; 
    } 
    else {
      kQfac3=1./3.; kQfac6=2./3.; 
      }  
    }
// **************** END IF BLOCK OF COLOR DECOMPOSITION ****************************************

  zeta1=cfac1*alphaQCD/vrel;
  zeta8=cfac8*alphaQCD/vrel;
  zeta3=cfac3*alphaQCD/vrel;
  zeta6=cfac6*alphaQCD/vrel;

  double SE = 1.; 
  
  if(kQfac1 > 0. || kQfac8 > 0. || kQfac3 > 0. || kQfac6 > 0.){ // Only if any of the prefactors is > 0, there is a nontrivial Sommerfeld factor. 
    SE=kQfac1*SommerfeldCoulomb(zeta1)+kQfac8*SommerfeldCoulomb(zeta8)+kQfac3*SommerfeldCoulomb(zeta3)+kQfac6*SommerfeldCoulomb(zeta6);
  }
  
  return SE;
}

void improveCrossSection(long n1,long n2,long n3,long n4 ,double PcmIn, double * res) {
  //model independent function for the Sommerfeld enhancement - see manual. 
  if(somm_flag==0){return;}  //No improvement of cross section

  else if(somm_flag==1){ //Sommerfeld enhancement for colored particles

// ****************** BEGIN s-wave Sommerfeld enhancement *******************  

/* The following code block decomposes the cross section of a process with
  two colored particles in the initial state into s- and p-wave parts. 
  It is model independent and NOT SUPPOSED TO BE CHANGED. 
  However, if the user wishes to use improveCrossSection in the original sense, he/she can
  use the "else" statement below to, for example, include loop improvements, as originally intended by
  the micrOMEGAs team.
  To use Sommerfeld enhancement AND loop improvements at the same time, use the additional 
  "else" statement and copy this code block there AND add your favorite loop improvements.  
*/
   
  // Initialization of informations regarding the particles participating in the process
  int nsub, nin, nout, err;
  int s1,s2,s3,s4;
  int q1,q2,q3,q4;
  int c1,c2,c3,c4;

  /* Do the names with proper pointers */
  char* name1 = new char[strlen(pdg2name(n1)) + 1]; // +1 for null terminator
  if (name1 != nullptr) {
    strcpy(name1, pdg2name(n1));
    // Use name1 as needed
  } 
  char* name2 = new char[strlen(pdg2name(n2)) + 1]; // +1 for null terminator
  if (name2 != nullptr) {
    strcpy(name2, pdg2name(n2));
    // Use name1 as needed
  } 
  char* name3 = new char[strlen(pdg2name(n3)) + 1]; // +1 for null terminator
  if (name3 != nullptr) {
    strcpy(name3, pdg2name(n3));
    // Use name1 as needed
  } 
  char* name4 = new char[strlen(pdg2name(n4)) + 1]; // +1 for null terminator
  if (name4 != nullptr) {
    strcpy(name4, pdg2name(n4));
    // Use name1 as needed
  } 

  qNumbers(name1, &s1,&q1,&c1);
  qNumbers(name2, &s2,&q2,&c2);
  qNumbers(name3, &s3,&q3,&c3);
  qNumbers(name4, &s4,&q4,&c4);

  if(c1!=1&&c2!=1){ //Enhancement only iff both particles are colored
    //Kinematics and couplings
    double mp1 = pMass(name1); double mp2 = pMass(name2); double mp3 = pMass(name3); double mp4 = pMass(name4); double mu = mp1*mp2/(mp1+mp2);
    double sqrts = (sqrt(pow(mp1,2)+pow(PcmIn,2))+sqrt(pow(mp2,2)+pow(PcmIn,2)));
    double vr= vRel(PcmIn,mp1,mp2); // relativistic
    double vcut=pow(10,-5); //to avoid the infinity at v=0
    double pcut=plowv(vcut,mp1,mp2); //corresponding p
    double QSE = mu*vr; // momentum of interaction
    double alphaSE = alphaQCD(QSE); // alpha strong evaluated at interaction momentum


    //Identifying the process in micrOMEGAs internally
    std::string procString = std::string(name1) + "," + name2 + "->" + name3 + "," + name4;
   	char* procName = new char [procString.length()+1];
    strcpy (procName, procString.c_str());
    numout* channel=newProcess(procName); //This gets not deallocated (MO internal reasons)
  	procInfo1(channel,&nsub,&nin,&nout);
    int nsub22_mem=nsub22; //Saving the matrix element and nsub22
    double (*sqme22_mem)(int, double, long double*, long double*, int*)=sqme22;

    //Sommerfeld enhancement
    double SE = SommerfeldFactor_BSMmodel(alphaSE,vr,c1,c2,c3,c4,s1,s2,s3,s4,n1,n2,n3,n4);

    // Calculation of the v=0 s-wave coefficient of sigma*vrel
    double swavecoeff =  cs22(channel,nsub,pcut,-1.0,+1.0,&err)/(3.89379E8)*vcut;

    double vcs = swavecoeff*SE;

    sqme22=sqme22_mem;
    nsub22=nsub22_mem;
  	*res=vcs/vr; //Dividing by vrel, since we improve only sigma, not sigma*v (vrel gets multiplied later automatically in the thermal average)

    delete[] procName;
  }
    
  delete[] name1;     
  delete[] name2;  
  delete[] name3;  
  delete[] name4;   
  return;
  }
  
  else{ // Implement your favorite loop improvements and/or Sommerfeld factors here
    return;
  }
}
