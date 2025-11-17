#include"SE_BSF_header.h"

int* getDarkSectorPDGs(int *num)
{
    /* This function returns an array of PDG numbers of ALL dark sector particles
       in the format: [A1, A2, ..., An, B1, -B1, B2, -B2, ..., Bm, -Bm] 
       where A_i are self-conjugate and B_j are non-self-conjugate */

    const int MAX_PARTICLES = 100;
    int tmpPDGs[2 * MAX_PARTICLES];  // Maximum possible size
    char txt[100];

    int err = sortOddParticles(txt);
    if (err) {
        printf("Error in sortOddParticles: %s\n", txt);
        *num = 0;
        return NULL;
    }

    int total_count = 0;
    int n = 0;
    double M;
    char *name;

    // First pass: collect self-conjugate particles
    //printf("Self-conjugate particles:\n");
    while ((name = nextOdd(n, &M)) && n < MAX_PARTICLES)
    {
        int spin2, charge3, cdim;
        int pdg = qNumbers(name, &spin2, &charge3, &cdim);
        if (pdg == 0) break;

        // Get antiparticle PDG to check if self-conjugate
        char *antiName = antiParticle(name);
        int spin2_anti, charge3_anti, cdim_anti;
        int pdg_anti = qNumbers(antiName, &spin2_anti, &charge3_anti, &cdim_anti);
        
        if (pdg == pdg_anti) {  // Self-conjugate
            tmpPDGs[total_count++] = pdg;
            //printf("  %s (PDG: %d, Mass: %.6f)\n", name, pdg, M);
        }
        n++;
    }

    // Second pass: collect non-self-conjugate particles as pairs
    //printf("\nNon-self-conjugate particle pairs:\n");
    n = 0;  // Reset counter
    while ((name = nextOdd(n, &M)) && n < MAX_PARTICLES)
    {
        int spin2, charge3, cdim;
        int pdg = qNumbers(name, &spin2, &charge3, &cdim);
        if (pdg == 0) break;

        // Get antiparticle PDG to check if non-self-conjugate
        char *antiName = antiParticle(name);
        int spin2_anti, charge3_anti, cdim_anti;
        int pdg_anti = qNumbers(antiName, &spin2_anti, &charge3_anti, &cdim_anti);
        
        if (pdg != pdg_anti && pdg > 0) {  // Non-self-conjugate and positive PDG (particle, not antiparticle)
            tmpPDGs[total_count++] = pdg;      // Particle
            tmpPDGs[total_count++] = pdg_anti; // Antiparticle
            //printf("  %s (PDG: %d) and %s (PDG: %d), Mass: %.6f\n", name, pdg, antiName, pdg_anti, M);
        }
        n++;
    }

    // Allocate exactly total_count entries
    int *pdgList = (int*)malloc(total_count * sizeof(int));
    if (!pdgList) {
        perror("malloc");
        *num = 0;
        return NULL;
    }

    for (int i = 0; i < total_count; i++)
        pdgList[i] = tmpPDGs[i];

    *num = total_count;
    return pdgList;
}

double smandel(double v, double m1, double m2){ //Mandelstam variable s of a given velocity v
  double s=pow(m1-m2,2)+2.*m1*m2*sqrt(1./(1.-pow(v,2)));
  return s;
}

double plowv(double v, double m1, double m2){ //relativistic correction to nonrelativistic relative velocity
  double p=m1*m2/sqrt((1.-pow(v,2))*(pow(m1,2)+pow(m2,2)+2.*m1*m2*sqrt(1./(1.-pow(v,2)))))*v;
  return p;
}

double smandelp(double Pcm, double m1, double m2){ //Mandelstam variable s of a given momentum p
  double  s=pow(sqrt(pow(m1,2)+pow(Pcm,2))+sqrt(pow(m2,2)+pow(Pcm,2)),2);
	return s;
}

double vRel(double Pcm, double m1, double m2){ //Relative velocity given a CM momentum Pcm
  double v=Pcm*sqrt(smandelp(Pcm,m1,m2))/(pow(Pcm,2)+sqrt((pow(m1,2)+pow(Pcm,2))*(pow(m2,2)+pow(Pcm,2))));
  return v;
}

double SommerfeldCoulomb(double zeta) {
  double s=1.;
  if (zeta!=0)  {   // Otherwise it has a singulartiy at the denominator and might go into troubles
    double h = 2.*M_PI*zeta;
    s = h/(1.-exp(-h));
    //s = exp(M_PI*zeta)*M_PI*zeta/(sinh(M_PI*zeta));
    //return s;
  }
  return s;
}

double BSF_XS_A(double T){ //Calculates the sum of perturbative annihilations and BSF cross sections 
  int Fast = 1;
  double Beps = 1.e-4;

  double pert_annihilation_cross_section = vSigmaA(T, Fast, Beps);
  double bound_state_contribution = 3.8937966E8*BoundStateCoannihilation(T); //Convert our result in Gev^{-2} to c*pb
  double total_annihilation_cross_section = pert_annihilation_cross_section + bound_state_contribution;

  return total_annihilation_cross_section;
}

double BSF_XS_S(double T){ //Calculates the perturbative semi-annihilations (no BSF decaying into dark sector are present)
  int Fast = 1;
  double Beps = 1.e-4;
  return vSigmaS(T, Fast, Beps);
}

double GammaDecayAnalytic(Parameters pars){ 
  /* Calculates the decay rate of an excited state according to eq. (62) of arXiv2112.01499v3. 
   */

  int s1,s2; //Quantum numbers of the coannihilators
  int q1,q2;
  int c1,c2;
  qNumbers(pdg2name(pars.n1), &s1,&q1,&c1);
  qNumbers(pdg2name(pars.n2), &s2,&q2,&c2);

  double spin_factor = 1.;
  if((s1 + s2) > 2){spin_factor = 0.5;}
  /* In the case of fermionic initial states, there is 
  an additional factor 1/2 in the <sigma*v>_s-wave. */

  double Gamma = 0.;
    
  int n_BSF = pars.n;

  //Definition of couplings
  double alpha_s_ann = pars.as_ann;
  double alpha_g_B = pars.ag_B;
  
  /* This block accounts for the thermal average */
  double mB = BS_mass(pars.m1, pars.m2, Ebin(pars.mu, alpha_g_B, n_BSF));
  double x_BS = mB/pars.T;
  double boost_factor = 1.; 
  double boost_test = bessK1(x_BS)/bessK2(x_BS);

  if(!isnan(boost_test)){boost_factor = boost_test;}

  if(pars.l==0){
    Gamma = spin_factor*boost_factor*pow(pars.mu*alpha_g_B/n_BSF, 3.)*alpha_s_ann*alpha_s_ann*casimir2(3)/(pars.m1*pars.m2); 
  }
  return Gamma;
}

double GammaDecayTriplet(Parameters pars){
  /* Here, the decay rate into 3 gluons for a spin triplet BS is calculates 
  according to eqs. (B9) of arXiv:2308.01336v2 and (67)  of arXiv:1702.01141v3 */
  double GammaSinglet = GammaDecayAnalytic(pars);
  double prefac = 4.*pars.as_ann*(M_PI*M_PI-9.)/(9.*M_PI);
  double Gamma = prefac*GammaSinglet;

  return Gamma;
}

double Ebin(double mu, double alpha, int n_BSF){ //Calculates the binding energy of a BS according to (48) of 2112.01499
  double alpha_eff = alpha/n_BSF;
  return 0.5*mu*alpha_eff*alpha_eff;
}

double BS_mass(double m1, double m2, double Ebin){ //Gives the mass of a bound state made of m1 and m2 bound by Ebin
  return m1 + m2 - Ebin;
}

double BS_dofs(int sb, int cb, int l){ //Returns the degrees of freedom of a BS with (micrOMEGA's) spin sb and color repr cb.
    int gB =  (2*l + 1.)*(sb + 1)*cb;
    return double(gB);
}

static void r_simpson(double(*func)(double, Parameters), Parameters pars, double *f, double a, double b, double eps, double *aEps, double *ans, double *aAns, int *deepness)
{
	double f1[5];
	int i;
	int d1 = *deepness + 1, d2 = *deepness + 1;
	double s1, s2, s3, e_err;

	s1 = (f[0] + 4 * f[4] + f[8])/6;
	s2 = (f[0] + 4 * f[2] + 2 * f[4] + 4 * f[6] + f[8]) / 12;
	s3 = (f[0] + 4 * f[1] + 2 * f[2] + 4 * f[3] + 2 * f[4] + 4 * f[5] + 2 * f[6] + 4 * f[7] + f[8]) / 24;

	if (!isfinite(s3))
	{
		*ans = s3;
		*aAns = s3;
		return;
	}

	e_err = eps * fabs(s3);
	i = 0;
	if ((fabs(s3 - s2) < e_err && fabs(s3 - s1) < 16 * e_err))
		i = 1;
	else if (fabs(s3 - s2) * (b - a) < 0.1 * (*aEps) && fabs(s3 - s1) * (b - a) < 1.6 * (*aEps))
	{
		i = 1;
		*aEps -= fabs((s3 - s2) * (b - a));
	}

	if (i || *deepness > 20)
	{
		*ans += s3 * (b - a);
		*aAns += (fabs(f[0]) + 4 * fabs(f[2]) + 2 * fabs(f[4]) + 4 * fabs(f[6]) + fabs(f[8])) * fabs(b - a) / 12;
		return;
	}

	for (i = 0; i < 5; i++)
		f1[i] = f[4 + i];
	for (i = 8; i > 0; i -= 2)
		f[i] = f[i / 2];
	for (i = 1; i < 8; i += 2)
		f[i] = (*func)(a + i * (b - a) / 16, pars);

	r_simpson(func, pars, f, a, (a + b) / 2, eps, aEps, ans, aAns, &d1);
	for (i = 0; i < 5; i++)
		f[2 * i] = f1[i];
	for (i = 1; i < 8; i+= 2)
		f[i] = (*func)((a + b) / 2 + i * (b - a) / 16, pars);
	r_simpson(func, pars, f, (a + b) / 2, b, eps, aEps, ans, aAns, &d2);
	if (d1 > d2)
		*deepness = d1;
	else
		*deepness = d2;
	return;
}

// Integration method from micrOMEGAs.
double simpsonArg(double (*func)(double, Parameters), Parameters pars, double a, double b, double eps)
{
	double f[9];
	double aEps;
	int i, j;

	aEps=0;
	if (a==b)
		return 0;
	for (i = 0; i < 9; i++)
	{
		f[i] = (*func) (a + i * (b - a) / 8, pars);
		aEps += fabs(f[i]);
	}
	if (aEps==0.)
		return 0;
	eps = eps / 2;
	aEps = eps * aEps * fabs(b - a) / 9;

	for (j = 0; ; j++)
	{
		double ans=0.0, aAns=0.0;
		int deepness = 1;
		r_simpson(func, pars, f, a, b, eps, &aEps, &ans, &aAns, &deepness);
		if (5 * aAns * eps > aEps || j >= 2 )
			return ans;
		if (!isfinite(aAns))
			return aAns;
		for (i = 0; i < 9; i++)
			f[i]=(*func)(a + i * (b - a) / 8, pars);
		aEps = aAns * eps;
	}
}

double Gamma_Ion_Milne_relation(double sigma_v_averaged, Parameters pars){ // Calculates <Gamma_ion> with the Milne relation according to (12) of 2112.01499
    /* Two particles n1 and n2 dissociate from a bound state with quantum numbers (n,l) and BSF rate sigma_v_averaged at temperature T. */
    double gdark1 = dofs_particle(pdg2name(pars.n1)); //Definition of parameters
    double gdark2 = dofs_particle(pdg2name(pars.n2));
    double alpha_g_B = pars.ag_B;
    double gB = BS_dofs(0, 1, pars.l); //ATTENTION! This assumes a spin 0 and color singlet bound state! Needs to be modified for the model at hand. 
    double binding_energy = Ebin(pars.mu, alpha_g_B, pars.n);
    double mB = BS_mass(pars.m1, pars.m2, binding_energy);
    
    return gdark1*gdark2/gB*pow(pars.T*pars.m1*pars.m2/(2.*M_PI*mB),1.5)*exp(-binding_energy/pars.T)*sigma_v_averaged;
}

double GammaIonGenMilne(double sigma_v_sum, Parameters pars){ 
  /* Calculates Gamma_Ion_Eff for a given sigma_v_sum. Derivation in my notes based on the definition in eq. (30) of arXiv:2112.01499v3 */
  
  int s1,s2; //Quantum numbers of the coannihilators
  int q1,q2;
  int c1,c2;
  qNumbers(pdg2name(pars.n1), &s1,&q1,&c1);
  qNumbers(pdg2name(pars.n2), &s2,&q2,&c2);

  double gdark1 = dofs_particle(pdg2name(pars.n1)); //Definition of parameters
  double gdark2 = dofs_particle(pdg2name(pars.n2));
  double alpha_g_B = pars.ag_B; //Value of alpha_g_B for the GS
  int nmax = pars.n;

  double Ebin_GS = Ebin(pars.mu, alpha_g_B, 1);
  double factor1 = gdark1*gdark2*pow(pars.mu*pars.T/(2.*M_PI), 1.5)*exp(-Ebin_GS/pars.T);
  double factor2 = 0.;
  double weighted_sum = 0.;

  for(int nBS = 1; nBS<=pars.n; nBS++){ //Loop over denominator quantities (see formula in my notes)
    double group_factor = pars.mu*alpha_group(c1, c2, 1, 1.); // ATTENTION! This is only correct for the singlet bound state configuration!!!
    double alpha_s_B = alphaQCD_binding(group_factor/double(nBS)); //earlier: alphaS_type(1,mu,0.1,alpha_table_331); // We choose a random vrel here since it is independent of it.
    double alpha_g_B_n = alpha_group(c1, c2, 1, alpha_s_B); 
    for(int lBS = 0; lBS<nBS; lBS++){
      double contrib = BS_dofs(0, 1, lBS)*exp((Ebin(pars.mu, alpha_g_B_n, nBS)-Ebin_GS)/pars.T);
      weighted_sum += contrib;
      if((s1 + s2) >= 2){weighted_sum += BS_dofs(2, 1, lBS)/BS_dofs(0, 1, lBS)*contrib;} //account for the spin triplets
    }
  }

  if(weighted_sum > 0){factor2 = 1./weighted_sum;}

  return factor1*factor2*sigma_v_sum;
}

double Sigma_v_Ionization_Equilibrium(Parameters pars){
  /* Calculates <sigma_BSF v_rel> according to eq. (33) of arXiv:2112.01499 */
  
  int s1,s2; //Quantum numbers of the coannihilators
  int q1,q2;
  int c1,c2;
  qNumbers(pdg2name(pars.n1), &s1,&q1,&c1);
  qNumbers(pdg2name(pars.n2), &s2,&q2,&c2);

  double gdark1 = dofs_particle(pdg2name(pars.n1)); //Definition of parameters
  double gdark2 = dofs_particle(pdg2name(pars.n2));
  double alpha_g_B = pars.ag_B;
  int nmax = pars.n;
  Parameters pars_loop = {pars.n1, pars.n2, 0, 0, pars.T, pars.as_ann, pars.ag_B, pars.m1, pars.m2, pars.mu}; //n, l and alpha_g_B will change in the loop
  
  double sum = 0.;
  double prefactor = pow(2.*M_PI/(pars.mu*pars.T), 1.5)/(gdark1*gdark2);
  
  for(int nBS = 1; nBS<=nmax; nBS++){
    double group_factor = pars.mu*alpha_group(c1, c2, 1, 1.); // ATTENTION! This is only correct for the singlet bound state configuration!!!
    double alpha_s_B = alphaQCD_binding(group_factor/double(nBS)); 
    double alpha_g_B_n = alpha_group(c1, c2, 1, alpha_s_B); 

    pars_loop.ag_B = alpha_g_B_n;

    double exponential_part = exp(Ebin(pars.mu, pars_loop.ag_B, nBS)/pars.T);
    pars_loop.n = nBS;
    for(int lBS = 0; lBS<nBS; lBS++){
      pars_loop.l = lBS;
      sum += BS_dofs(0, 1, lBS)*exponential_part*GammaDecayAnalytic(pars_loop);
      if((s1 + s2) >= 2){sum += BS_dofs(2, 1, lBS)*exponential_part*GammaDecayTriplet(pars_loop);} //Account for the triplet part
    }
  }
  return prefactor*sum;
}

double Gamma_Decay_effective_ET(Parameters pars){ 
  /* Calculate the effective decay rate in the Efficient Transition limit based on 
  eq. (12) of arXiv:2112.01499 and simplified in my notes*/
  
  int s1,s2; //Quantum numbers of the coannihilators
  int q1,q2;
  int c1,c2;
  qNumbers(pdg2name(pars.n1), &s1,&q1,&c1);
  qNumbers(pdg2name(pars.n1), &s2,&q2,&c2);

  double gdark1 = dofs_particle(pdg2name(pars.n1)); //Definition of parameters
  double gdark2 = dofs_particle(pdg2name(pars.n2));
  double alpha_g_B = pars.ag_B; //Value for alpha_g_B for the ground state
  int nmax = pars.n;

  double Ebin_GS = Ebin(pars.mu, alpha_g_B, 1); //Values for the ground state
  double mB_GS = BS_mass(pars.m1, pars.m2, Ebin_GS);
  double numerator = 0.;
  double denominator = 0.;

  double result = 0.; //Useful parameters
  double dof_i = 0.;
  double dof_i_triplet = 0.;
  double mBi = 0.;
  double mBi_ratio = 0.;
  double Ebi = 0.;
  double exponential_factor = 0.;
  double common_factor = 0.;
  double common_factor_triplet = 0.;
  Parameters pars_loop = {pars.n1, pars.n2, 0, 0, pars.T, pars.as_ann, pars.ag_B, pars.m1, pars.m2, pars.mu}; //n, l and alpha_g_B will change in the loop

  for(int nBS = 1; nBS<=nmax; nBS++){ //Masterloop 
    double group_factor = pars.mu*alpha_group(c1, c2, 1, 1.); // ATTENTION! This is only correct for the singlet bound state configuration!!!
    double alpha_s_B = alphaQCD_binding(group_factor/double(nBS)); 
    double alpha_g_B_n = alpha_group(c1, c2, 1, alpha_s_B); 

    pars_loop.ag_B = alpha_g_B_n;
    Ebi = Ebin(pars.mu,alpha_g_B_n, nBS);
    mBi = BS_mass(pars.m1, pars.m2, Ebi);
    mBi_ratio = pow(mBi/mB_GS, 1.5);
    exponential_factor = exp((mB_GS - mBi)/pars.T);
    pars_loop.n = nBS;
    for(int lBS = 0; lBS<nBS; lBS++){
      dof_i = BS_dofs(0, 1, lBS);
      if((s1 + s2) >= 2){dof_i_triplet = BS_dofs(2, 1, lBS);}
      pars_loop.l = lBS;
      common_factor = dof_i*mBi_ratio*exponential_factor;
      common_factor_triplet = dof_i_triplet*mBi_ratio*exponential_factor;
      denominator += (common_factor + common_factor_triplet);
      numerator += (common_factor*GammaDecayAnalytic(pars_loop) + common_factor_triplet*GammaDecayTriplet(pars_loop));
    }
  }

  if(denominator>0.){result = numerator/denominator;}

  return result;
}

double Pnl2_prefac_squared(double zeta_s, double zeta_bt, int n, int l){ // residual irreducible factor of (P(n, l))^2 from eq. (A6) of arXiv:2308.01336v1 
    double prefac = 1/double(n);
    double prod = 1.;
    double ratio = 4./(1./zeta_bt + zeta_bt);
    double ratio2 = ratio*ratio;
    double n2 = double(n*n);
    double zeta_s2 = zeta_s*zeta_s;
    if(l>0){
        for(int j = 0; j < l; j++){
            double factor = ratio2*(double(j*j) + zeta_s2)/(n2 - double((j+1)*(j+1)));
            prod = prod*factor;
        }
    }
    
    return prefac*prod;
}

double IR_bracket_part_high(double zeta_s, double zeta_bt, std::array<double, 5> Rnl, int l){ //SQUARED bracket part of eq. (A4) of arXiv:2308.01336v1 for xmax = n-l.
    /* The exp(atan(zeta_bt)) factor is here for converginience (== convergence && convenience)! */
    double prefac = exp(2.*zeta_s*atan(zeta_bt));
    double h1 = double(l + 2); //Auxiliary quantities
    double h2 = 2.*zeta_s;
    double bracket_part = (h1*zeta_bt + zeta_s)*Rnl[0] - 2.*(h1*zeta_bt + h2)*Rnl[1] + 6.*zeta_s*Rnl[2] + 2.*(h1*zeta_bt - h2)*Rnl[3] + (zeta_s - h1*zeta_bt)*Rnl[4];
    double product = prefac*abs(bracket_part);
    return product*product;
}

double IR_bracket_part_low(double zeta_s, double zeta_bt, std::array<double, 5> Rnl, int l){ //SQUARED bracket part of eq. (A5) of arXiv:2308.01336v1 for xmax = n-l.
    /* The exp(atan(zeta_bt)) factor is here for converginience (== convergence && convenience)! */
    double prefac = exp(2.*zeta_s*atan(zeta_bt));
    double ld = double(l);
    double h1 = 1. + ld; //Auxiliary quantities
    double h2 = ld + h1;
    double zs2 = zeta_s*zeta_s;
    double zs3 = zeta_s*zs2;
    double zb2 = zeta_bt*zeta_bt;
    double part0, part1, part2, part3, part4;
    part0 = (ld*h1*zeta_bt*(-3. + h2)*zb2) + (-1. - 3.*ld + 3.*h1*h2*zb2)*zeta_s + 6.*h1*zeta_bt*zs2 + 2.*zs3;
    part1 = ld*h1*zeta_bt*(3. + h2*zb2) + 2.*zeta_s + 6.*ld*zeta_s - 6.*h1*zeta_bt*zs2 - 4.*zs3; 
    part2 = 6.*zeta_s*(1. + 3.*ld + h1*h2*zb2 - 2.*zs2);
    part3 = ld*h1*zeta_bt*(3. + h2*zb2) - 2.*(ld + h2)*zeta_s - 6.*h1*zeta_bt*zs2 + 4.*zs3;
    part4 = ld*h1*zeta_bt*(-3. + h2*zb2) + zeta_s + 3.*ld*zeta_s - 3.*h1*h2*zb2*zeta_s + 6.*h1*zeta_bt*zs2 - 2.*zs3;
    double bracket_part = part4*Rnl[4] + 2.*part3*Rnl[3] + part2*Rnl[2] - 2.*part1*Rnl[1] - part0*Rnl[0];
    double product = prefac*abs(bracket_part);
    return product*product;
}

double bsf_cs_integrand_excited_states(double vrel, Parameters pars){
    double bsf_integrand = 0.;
    double bsf_sigmav = 0.;
    
    int s1,s2; //Quantum numbers of the coannihilators
    int q1,q2;
    int c1,c2;
    qNumbers(pdg2name(pars.n1), &s1,&q1,&c1);
    qNumbers(pdg2name(pars.n2), &s2,&q2,&c2);
    
    int n_BSF = pars.n;
    int l_BSF = pars.l;
    double temp = pars.T;
    
    double m1 = pars.m1; //Masses of the two coannihilators
    double m2 = pars.m2;
    double mu = pars.mu; //reduced mass of the bound state
    double mu_over_2T = mu/(2.*temp);
    double gdark1 = dofs_particle(pdg2name(pars.n1)); //dofs of the two coannihilators
    double gdark2 = dofs_particle(pdg2name(pars.n2));

    //Definition of couplings
    double alpha_s_ann =  pars.as_ann;
    double alpha_g_B = pars.ag_B;
    double alpha_g_S = alpha_group(c1, c2, 8, alpha_s_ann); // ATTENTION! we have to use the octet scatteriung state as we evaluate 8->1
    
    //Important factors
    double cas2 = casimir2(c1); // = CF = 4/3
    double dR = double(c1); // = Nc = 3
    double zeta_s = alpha_g_S/vrel;
    double zeta_b = alpha_g_B/vrel;
    double alpha_B_eff = double(alpha_g_B/n_BSF);
    double zeta_b_tilde = double(zeta_b/n_BSF);
    double omega = mu/2.*(vrel*vrel + alpha_B_eff*alpha_B_eff); 
    double alpha_s_BSF = alphaQCD(omega); 
    double bose = 1./(1. - exp(-omega/temp));
    double boltzmann = pow(mu_over_2T/M_PI, 3./2.)*exp(-mu_over_2T*pow(vrel,2.));
    double intmeasure = 4.*M_PI*vrel; // we extract one power of v to cancel the 1/v from SBSF for the computer for v->0
    double S_BSFv = 0.; 
    
    double ld = double(l_BSF); //for the upcoming numerics
    double high_l_contribution = 0.;
    double low_l_contribution = 0.;
    //Auxiliary factors
    double h1 = 1. + ld;
    double zs2 = zeta_s*zeta_s;
    
    if(vrel>0.){ // to avoid the 0/0 for v=0. In the limit vrel -> 0, sigmaBSF*vrel -> 0.
        /* Here, eq. (A1) of arXiv:2308.01336v1 is performed  */
        std::array<double, 5> Rnl_array = Rnl_recursion(n_BSF, l_BSF, zeta_b_tilde, zeta_s);
        high_l_contribution = 4.*(ld*ld + zs2)*h1*(zs2 + h1*h1)*IR_bracket_part_high(zeta_s, zeta_b_tilde, Rnl_array, l_BSF);
        if(l_BSF>0){
            low_l_contribution = ld*IR_bracket_part_low(zeta_s, zeta_b_tilde, Rnl_array, l_BSF);
        }
        S_BSFv = Pnl2_prefac_squared(zeta_s, zeta_b_tilde, n_BSF, l_BSF)*(high_l_contribution + low_l_contribution);
        double prefac = 64.*M_PI*M_PI*cas2*alpha_s_BSF*alpha_g_B*alpha_g_B/(3.*dR*dR*mu*mu)*(zeta_b_tilde/zeta_s)/(exp(2.*M_PI*zeta_s)-1.)/pow(double(n_BSF)*(1. + zeta_b_tilde*zeta_b_tilde), 3.);
        bsf_sigmav = prefac*S_BSFv;
    }
        
    bsf_integrand = intmeasure*boltzmann*bose*bsf_sigmav;
    return bsf_integrand;
}

std::array<double, 5> Rnl_recursion(int n, int l, double zeta_bt, double zeta_s){
  /* Returns an array of values {Rnl(n-l-i)}, i = 1,...,5 according to eq. (A7) of arXiv:2308.01336v2. */
    std::array<double, 5> R_array; //Final output array containing only the 5 needed elements 
    int xmax = n - l + 4; //number of entries to be calculated in the large array
    int xrelevant = xmax - 5; //index at which the relevant entries of the large array start 
    double R_array_long[xmax]; //Initialization of the full array containing ALL elements
    
    //Definition of constants
    double zb2 = zeta_bt*zeta_bt;
    double h1 = zb2 - 1.;
    double h2 = zb2 + 1.;
    double h3 = 4.*zeta_bt*zeta_s;
    double ld = double(l);

    //Recursion steps
    for(int i = 0; i<xmax; i++){R_array_long[i] = 0;} //Initial values 0.
    R_array_long[4] = 1.;
    R_array_long[5] = (h3 - 2.*(3. + ld)*h1)/h2;
    for(int j = 6; j<xmax; j++){
      double xrec = double(j) - 4.;
      R_array_long[j] = (h3 - 2.*(2. + ld + xrec)*h1)/(xrec*h2)*R_array_long[j-1] - (4. + 2.*ld + xrec)/xrec*R_array_long[j-2];
    }

    for(int i = 0; i<5; i++){ //Final output array declaration
      R_array[i] = R_array_long[xrelevant + i];      
    }

    return R_array; //array returned
}

double casimir2(int color) // = CF = (Nc^2-1)/(2*Nc); CalcHEP up to 3.8.7 implements only 1, 3, 6 and 8 representations
{
    switch (color)
    {
      case 1: return 0.;
      case 3: return 4./3.;
      case -3: return 4./3.;
		  case 6: return 10./3.;
		  case 8: return 3;
  //  case 10: return 6;
  //  case 27: return 8;
		default:
			printf("WARNING: casimir2 called with color: %d.\n", color);
			exit(42);
    }
}

double alpha_group(int r1, int r2, int r, double alpha_s){ 
    /* Returns the coupling alpha_g for a given alpha_s according to eq.(2.14) of arXiv:1805.01200v2 */
    return 0.5*(casimir2(r1) + casimir2(r2) -casimir2(r))*alpha_s;        
}

double alphaQCD_binding(double fac){
    /* This function returns alpha_s^B as defined in table 1 of arXiv:1805.01200v2. It uses an iterative method.  */
    
    double diff = 10.;
    double alpha_old = 1.; //initial guess for the correct scale.
    double alpha_new, alpha_s_B = 0.;
    
    while(diff > 0.001){
      alpha_new = alphaQCD(fac*alpha_old);
      diff = abs(alpha_new - alpha_old);
      alpha_old = alpha_new; 
    }
    alpha_s_B = alpha_new;
    
    return alpha_s_B;
}

double BSF_XS_TA_ex(long n1, long n2, double T, int scenario, int n_BSF){ 
  /* Thermal average of the full BSF contribution between two particles n1 and n2. 
  T is the temperature, scenario is the BSF scenario (no/efficient transition limit, ionization equilibrium, full or no solution) and 
  n_BSF is the number of included states. */
  if((scenario==0)||(n_BSF==0)){return 0.;}
  int s1,s2; //Quantum numbers of the coannihilators
  int q1,q2;
  int c1,c2;
  qNumbers(pdg2name(n1), &s1,&q1,&c1);
  qNumbers(pdg2name(n2), &s2,&q2,&c2);

  double total_spin = double((s1 + 1)*(s2 + 1)); //Total spin normalization
  bool triplet = total_spin >= 2;

  double m1 = pMass(pdg2name(n1)); //Masses of the two coannihilators
  double m2 = pMass(pdg2name(n2));
  double mu = m1*m2/(m1+m2); //reduced mass of the bound state

  //Definition of couplings
  //double ann_scale = (2*Mcdm + T)/3; //Value used for GGScale in micrOMEGAs
  double ann_scale = sqrt(1.5*mu*T); //Scale relevant for annihilation processes
  double alpha_s_ann = alphaQCD(ann_scale);
  // actually we would have to evaluate it at Q = mu*vrel.
  // We use the averaged velocity vrel = sqrt(3T/MY) -> Q = sqrt(3 MY T)
  
  /* Calculation of alpha_s_B for every n */
  double alpha_g_B_vec[n_BSF];
  double alpha_s_B_group_factor = mu*alpha_group(c1, c2, 1, 1.); // ATTENTION! This is only correct for the singlet bound state configuration!!!
  for(int i = 0; i<n_BSF; i++){
    double alpha_s_B = alphaQCD_binding(alpha_s_B_group_factor/double(i+1)); 
    double alpha_g_B = alpha_group(c1, c2, 1, alpha_s_B); 
    alpha_g_B_vec[i] = alpha_g_B;
  }
  
  double bsf_eff = 0.; //Initialize final output of BSF cross section. 
  int total_number_of_states = (n_BSF*(n_BSF + 1))/2; //l = 0, 1, ..., n-1
  double bsfaverages[total_number_of_states]; // final array of <sigma v> for every bound state (spin singlet)
  double Ion_rates[total_number_of_states]; //array of ionisation rates for every bound state
  double Decay_rates[total_number_of_states];// Array of decay_rates 

  double bsfaverages_triplet[total_number_of_states]; // final array of <sigma v> for every bound state (spin triplet)
  double Ion_rates_triplet[total_number_of_states]; //array of ionisation rates for every bound state (spin triplet)
  double Decay_rates_triplet[total_number_of_states];// Array of decay_rates (spin triplet)

  for(int i = 0; i<total_number_of_states; i++){//Initializations
    bsfaverages[i] = 0.;
    Ion_rates[i] = 0.;
    Decay_rates[i] = 0;
    bsfaverages_triplet[i] = 0.;
    Ion_rates_triplet[i] = 0.;
    Decay_rates_triplet[i] = 0.; 
  }

  switch(scenario) { 
    case 1:   //No transition limit
      if((c1==3 && c2==-3) || (c1==-3 && c2==3)){ //Only fundamental-antifundamental BS up to now 
        
        int master_index = 0;
        for(int n = 1; n<=n_BSF; n++){//Master loop: Calculate ALL rates
          double alpha_g_B_n = alpha_g_B_vec[n-1]; 
          for(int l = 0; l<n; l++){
            Parameters pars = {n1, n2, n, l, T, alpha_s_ann, alpha_g_B_n, m1, m2, mu}; //Parameters relevant for the calculation of <sigma*v>, GammaIon and GammaDecay.
            Decay_rates[master_index] = GammaDecayAnalytic(pars);
            if(Decay_rates[master_index] > 0.){
              bsfaverages[master_index] = simpsonArg(bsf_cs_integrand_excited_states, pars, 0, 1., 0.00001)/total_spin;
              Ion_rates[master_index] = Gamma_Ion_Milne_relation(bsfaverages[master_index], pars);
            }

            if(triplet){
              Decay_rates_triplet[master_index] = GammaDecayTriplet(pars);
              if(Decay_rates_triplet[master_index] > 0.){
                bsfaverages_triplet[master_index] = 3.*bsfaverages[master_index]; //factor 3 coming from sum over final states.
                Ion_rates_triplet[master_index] = Ion_rates[master_index]; //prefactor coming from different BS dofs cancels with multiplicity of sigma*v. 
              }
            }

            master_index += 1;
          }
        }
        
        for(int i = 0; i<total_number_of_states; i++){ //Final step of adding the contributions together
            if(Decay_rates[i] > 0.){
              bsf_eff += bsfaverages[i]/(1. + Ion_rates[i]/Decay_rates[i]);
            }
            if(Decay_rates_triplet[i] > 0.){
              bsf_eff += bsfaverages_triplet[i]/(1. + Ion_rates_triplet[i]/Decay_rates_triplet[i]);
            }
        }
      }
      else{
        bsf_eff = 0.;
      }
      break;

    case 2:  //Efficient transition limit
      if((c1==3 && c2==-3) || (c1==-3 && c2==3)){ //Only fundamental-antifundamental BS up to now 
        double bsf_sum = 0.; //Definition of final ingredients
        double GammaDecayEff = 0.;
        double GammaIonEff = 0.;

        for(int n = 1; n<=n_BSF; n++){ //Calculation of <sigmaBSF_vrel>_sum
          for(int l = 0; l<n; l++){
            Parameters pars = {n1, n2, n, l, T, alpha_s_ann, alpha_g_B_vec[n-1], m1, m2, mu}; //Parameters relevant for the calculation of <sigma*v>, GammaIon and GammaDecay.
            double bsf_contrib = simpsonArg(bsf_cs_integrand_excited_states, pars, 0, 1., 0.00001)/total_spin;
            bsf_sum += bsf_contrib;
            if(triplet){bsf_sum += 3.*bsf_contrib;} //add the triplet contribution    
          }
        }

        //Calculation of Gamma_Ion/Dec_Eff
        Parameters pars_global = {n1, n2, n_BSF, 0, T, alpha_s_ann, alpha_g_B_vec[0], m1, m2, mu}; 
        /* The value of l is irrelevant here, since we sum over it in the internal function. 
        alpha_g_B is for the GS, will be modified inside the functions below accordingly. */
        GammaIonEff = GammaIonGenMilne(bsf_sum, pars_global);
        GammaDecayEff = Gamma_Decay_effective_ET(pars_global);

        if(GammaDecayEff>0){ //Final step according to eq. (31) of arXiv:2112.01499v3
          bsf_eff = bsf_sum/(1. + GammaIonEff/GammaDecayEff);
        }
      }
      else{
        bsf_eff = 0.;
      }
      break;

    case 3:   //Ionization equilibrium
      if((c1==3 && c2==-3) || (c1==-3 && c2==3)){ //Only fundamental-antifundamental BS up to now 
        Parameters pars_global = {n1, n2, n_BSF, 0, T, alpha_s_ann, alpha_g_B_vec[0], m1, m2, mu}; 
        /* The value of l is irrelevant here, since we sum over it in the internal function 
        alpha_g_B is for the GS and is adapted accordingly inside the function below. */
        bsf_eff = Sigma_v_Ionization_Equilibrium(pars_global);
      }
      else{
        bsf_eff = 0.;
      }
      break;
    
    case 4:  //Full solution
      if((c1==3 && c2==-3) || (c1==-3 && c2==3)){ //Only fundamental-antifundamental BS up to now
        
        char* QEDCouplingName[5] = {nullptr}; // Initialize all elements to nullptr
        // Allocate and initialize CouplingName
        const char* tempName = "aEWM1";
        QEDCouplingName[0] = (char*)malloc((strlen(tempName) + 1) * sizeof(char));
        if (QEDCouplingName[0] == nullptr) {
          // Handle allocation failure
          return 1;
        }
        strcpy(QEDCouplingName[0], tempName);

        double alphaQED = 1./findValW(QEDCouplingName[0]); //We do not include a running QED coupling (yet). 
        free(QEDCouplingName[0]); //Free the allocated memory. 

        double Gamma_BS_total[total_number_of_states]; // GammaIon_i + GammaDec_i + sum_{j \neq i} GammaTrans_ji as defined in eq. (19) in arXiv:2112.01499
        double Gamma_array[total_number_of_states]; // GammaIon_j/Gamma_BS_total_j
        const int total_num_trans = total_number_of_states*total_number_of_states; //Size of transition matrix
        double Gamma_trans_matrix[total_num_trans]; //Matrix of transitions
        double MMatrix[total_num_trans]; //Matrix M defined in eq. (21) in arXiv:2112.01499
        double InverseMMatrix[total_num_trans]; //Inverse of Matrix M

        /* Initialize respective quantities for the triplet states */
        double Gamma_BS_total_triplet[total_number_of_states]; // GammaIon_i + GammaDec_i + sum_{j \neq i} GammaTrans_ji as defined in eq. (19) in arXiv:2112.01499
        double Gamma_array_triplet[total_number_of_states]; // GammaIon_j/Gamma_BS_total_j
        double MMatrix_triplet[total_num_trans]; //Matrix M defined in eq. (21) in arXiv:2112.01499
        double InverseMMatrix_triplet[total_num_trans]; //Inverse of Matrix M

        int master_index = 0;
        for(int n = 1; n<=n_BSF; n++){//Master loop: Calculate ALL rates and initialize other quantities
          for(int l = 0; l<n; l++){
            Parameters pars = {n1, n2, n, l, T, alpha_s_ann, alpha_g_B_vec[n-1], m1, m2, mu}; //Parameters relevant for the calculation of <sigma*v>, GammaIon and GammaDecay.
            Decay_rates[master_index] = GammaDecayAnalytic(pars);
            double sigmav_temp = simpsonArg(bsf_cs_integrand_excited_states, pars, 0, 1., 0.00001)/total_spin;
            double GammaIon_temp = Gamma_Ion_Milne_relation(sigmav_temp, pars);
            bsfaverages[master_index] = sigmav_temp;
            Ion_rates[master_index] = GammaIon_temp;
            Gamma_BS_total[master_index] = 0.; //Initializations
            Gamma_array[master_index] = 0.;

            //Initialize/calculate respective rates for the triplet
            Gamma_BS_total_triplet[master_index] = 0.;
            Gamma_array_triplet[master_index] = 0.;
            
            if(triplet){ 
              Decay_rates_triplet[master_index] = GammaDecayTriplet(pars);
              bsfaverages_triplet[master_index] = 3.*sigmav_temp;
              Ion_rates_triplet[master_index] = GammaIon_temp;
            }
            
            master_index += 1;
          }
        }
        
        // Define the transition rates here 
        for(int i = 0; i<total_num_trans; i++){//Initializations
          Gamma_trans_matrix[i] = 0.;
          MMatrix[i] = 1.;
          InverseMMatrix[i] = 0.;
          MMatrix_triplet[i] = 1.;
          InverseMMatrix_triplet[i] = 0.;
          }
        
        /* Calculating gamma trans matrix elements */
        double trans_prefac = 4.*abs(double(q1)*double(q2))*alphaQED/27.; //constant prefactor of all transition rates 
        for(int i = total_number_of_states-1; i>=0; i--){
          for(int j = i-1; j>=0; j--){
            std::array<int, 2> nl_i = index_conversion(i); //Quantum numbers of (higher n') i state state and (lower n) j state
            std::array<int, 2> nl_j = index_conversion(j);

            if((nl_i[0]!=nl_j[0]) && abs(nl_i[1] - nl_j[1])==1){ //Only calculate matrix element if n != n' and \Delta l = \pm 1 (assuming Enl = Enl' \forall l, l')
              double alpha_i = alpha_g_B_vec[nl_i[0]-1]; //alpha_g_B for BS i. 
              double alpha_j = alpha_g_B_vec[nl_j[0]-1]; //alpha_g_B for BS j. 
              //Binding energies and masses of (higher) state i and (lower) state j
              double EBi = Ebin(mu, alpha_i, nl_i[0]); 
              double EBj = Ebin(mu, alpha_j, nl_j[0]);         
              double mBi = BS_mass(m1, m2, EBi);
              double mBj = BS_mass(m1, m2, EBj);
              double omega = EBi - EBj;
              
              TransitionParameters parsT = {nl_j[0], nl_j[1], nl_i[0], nl_i[1], mu, alpha_j, alpha_i};
              double Gamma_trans_ij = trans_prefac*(2.*double(nl_j[1]) + 1.)*pow(abs(omega), 3.)*GammaTrans(parsT); //Transition rate from (higher n') state i to (lower n) state j  
              /* Ratio of Gamma_ij/Gamma_ji according to the Milne relation (13) in arXiv:2112.01499v2 */
              double Milne_factor = pow(mBi/mBj, 1.5)*exp(omega/T); 
              
              Gamma_trans_matrix[i*total_number_of_states + j] = Gamma_trans_ij;
              Gamma_trans_matrix[j*total_number_of_states + i] = Milne_factor*Gamma_trans_ij;
            }
          }
        }
        
        // Calculate Gamma_BS_total and Gamma_array 
        for(int i = 0; i<total_number_of_states; i++){
          double trans_addendum = 0.;
          for(int j = 0; j<total_number_of_states; j++){
            trans_addendum += Gamma_trans_matrix[total_number_of_states*i + j]; //Sum over rows: sum_j GammaTrans_{i->j}
          }
          double sum = Ion_rates[i] + Decay_rates[i] + trans_addendum;
          double sum_triplet = Ion_rates_triplet[i] + Decay_rates_triplet[i] + trans_addendum;

          if(sum>0.){
          Gamma_BS_total[i] = sum;
          Gamma_array[i] = Ion_rates[i]/sum;
          }

          if(triplet && (sum_triplet > 0.)){
            Gamma_BS_total_triplet[i] = sum_triplet;
            Gamma_array_triplet[i] = Ion_rates_triplet[i]/sum_triplet;
          }
        }

        // Set up Matrix M according to eq. (21) in arXiv:2112.01499 
        for(int i = 0; i<total_number_of_states; i++){
          for(int j = i+1; j<total_number_of_states; j++){
            MMatrix[i*total_number_of_states + j] = -Gamma_trans_matrix[i*total_number_of_states + j]/Gamma_BS_total[i];
            MMatrix[j*total_number_of_states + i] = -Gamma_trans_matrix[j*total_number_of_states + i]/Gamma_BS_total[j];
            
            if(triplet){
              MMatrix_triplet[i*total_number_of_states + j] = -Gamma_trans_matrix[i*total_number_of_states + j]/Gamma_BS_total_triplet[i];
              MMatrix_triplet[j*total_number_of_states + i] = -Gamma_trans_matrix[j*total_number_of_states + i]/Gamma_BS_total_triplet[j];
            }
            
          }
        }
        // MATRIX INVERSION 
        /* We need to convert MMatrix into long double */
         long double *Mat = (long double*)malloc(total_num_trans * sizeof(long double));
    	  if (Mat == NULL) {
          fprintf(stderr, "Memory allocation failed for Mat\n");
          return 1;
         }

        for(int i = 0; i<total_num_trans; i++){
          Mat[i] = double(MMatrix[i]);
        }
        long double *Mij = Mat;
        double* MMatrixInverse = MatrixInversion(total_number_of_states, Mij);

        /* Very final step of adding the contributions together according to eq. (23) in arXiv:2112.01499  */
        for(int i = 0; i<total_number_of_states; i++){ 
          double factor = 0.;
          for(int j = 0; j<total_number_of_states; j++){
            factor += MMatrixInverse[total_number_of_states*i + j]*Gamma_array[j];          
          }
          bsf_eff += bsfaverages[i]*(1. - factor);         
        }

        if(triplet){//Redo the steps for the triplet
          long double *Mat_triplet = (long double*)malloc(total_num_trans * sizeof(long double));
    	    if (Mat_triplet == NULL) {
            fprintf(stderr, "Memory allocation failed for Mat\n");
            return 1;
          }
          for(int i = 0; i<total_num_trans; i++){
            Mat_triplet[i] = double(MMatrix_triplet[i]);
          }
        long double *Mij_triplet = Mat_triplet;
        double* MMatrixInverse_triplet = MatrixInversion(total_number_of_states, Mij_triplet);

        /* Very final step of adding the contributions together according to eq. (23) in arXiv:2112.01499  */
        for(int i = 0; i<total_number_of_states; i++){
            double factor_triplet = 0.;
            for(int j = 0; j<total_number_of_states; j++){
              factor_triplet += MMatrixInverse_triplet[total_number_of_states*i + j]*Gamma_array_triplet[j];
            }
            bsf_eff += bsfaverages_triplet[i]*(1. - factor_triplet);
          }
          if (MMatrixInverse_triplet != NULL) {
          free(MMatrixInverse_triplet);
          }
          
        free(Mat_triplet); 
        }
        
        if (MMatrixInverse != NULL) {
          free(MMatrixInverse);
          }
        
        free(Mat); 
      }
      else{
        bsf_eff = 0.;
      }
      break;
    
    default:
      std::cout << "No proper BSF scenario chosen" << std::endl;
  }
  if(isnan(bsf_eff)){std::cout << "bsf_eff = NAN: T = " << T << std::endl;}
  return bsf_eff;
}

double GammaTrans(TransitionParameters pars){
    /* Calculates the expression of eq. (A8) of arXiv:2308.01336v2 */
	double Gamma = 0.; //Initialize final output
	
	double nd = double(pars.n); //Convert quantum numbers to double
	double npd = double(pars.np);
	double ld = double(pars.l);
	double lpd = double(pars.lp);
	
	double kappa = pars.alpha/nd; //Initialize important quantities (we factor out mu)
	double kappa_p = pars.alpha_p/npd;
	
	/* Define many auxiliary quantities */
	double h1 = kappa + kappa_p;
	double h2 = h1*h1;
	double z = 4.*kappa*kappa_p/h2;
	double h3 = ld + 1.; 
	double h4 = kappa - kappa_p;
	double h5 = h3 - nd; //l + 1 -n
	double h6 = 2.*h3; //2*l + 2
	double h7 = h4*h4; //(k - k')^2
	double h8 = pars.alpha - pars.alpha_p; 
  double h9 = h3 + npd; // n'+l+1
  double h10 = nd + ld;
  double h11 = 2.*h8/(h1*h4); //2*(alpha - alpha')/(k^2 - k'^2)
  double h12 = h6 - 2.; //2*l
  double h13 = h9-1.; //n' + l
  double h14 = npd - nd;
  double h15 = pow(1.-z, h14); //(1-z)^(n' - n)
  double h16 = h13 - 1.; //n' + l - 1
  double h17 = h6 - 1.; // 2*l+1
  double h18 = pars.alpha/pars.alpha_p;
	
	double prefac = 0.; //Initialize quantities important for the final output
	double prod1 = 1.; 
	double prod2 = 1.;
	double bracket_part = 0.;
    
    switch(pars.l - pars.lp){
		case -1: // l' = l + 1
			/* Implementation of the square of eq.(A12) of arXiv:2308.01336v2 (derivation in my notes) */
			prefac = h18*h3*h15*pow(z, h6)*(npd + h3)*(npd - h3)*(npd - ld)*(nd + ld)/(nd*nd*h17*h17);
            
      for(int k = 1; k<=2*pars.l; k++){ //Calculate prod1 and prod2
          double kd = double(k);
          prod1 = prod1*(h9/kd - 1.);
          prod2 = prod2*(h10/kd - 1.);                
      }
			
			bracket_part += (kappa_p - h8)/h2*hypergeometric_prior(h5, h9 + 1, h6, z); 
			bracket_part += h11*hypergeometric_prior(h5, h9, h6, z); //The prefactor of this copntribution is zero
			bracket_part -= (kappa_p + h8)/h7*hypergeometric_prior(h5, h13, h6, z); 
			break;

		case 1: // l' = l - 1
      /* Implementation of the square of eq.(A13) of arXiv:2308.01336v2 (derivation in my notes) */
			prefac = ld*h15*pow(z, h12)*h10*(nd-ld)*h16/(h18*npd*npd*(npd - ld));
            
      for(int k = 1; k<2*pars.l; k++){ //Calculate prod1 and prod2
          double kd = double(k);
          prod1 = prod1*(h10/kd - 1.);
          prod2 = prod2*(h16/kd - 1.);                
      }
			
			bracket_part -= (kappa + h8)/h7*hypergeometric_prior(h5 - 2., h10, h12, z); 
			bracket_part -= h11*hypergeometric_prior(h5-1., h9-1., h12, z); //The prefactor of this copntribution is zero
			bracket_part += (kappa - h8)/h2*hypergeometric_prior(h5, h13, h12, z); 
			break;

		default: 
			std::cout << "Selection rule violated, Gamma_trans = 0" << std::endl;
			break;
	}
	
  bracket_part = bracket_part/pars.mu;
	Gamma = prefac*prod1*prod2*bracket_part*bracket_part;
  if(Gamma<0.){std::cout << "GammaTrans is negative! check function GammaTrans again!" << std::endl;}
  if(isnan(Gamma)){std::cout << "GammaTrans = " <<  Gamma << std::endl;}
	return Gamma;
}

std::array<int, 2> index_conversion(int i){
  /* Converts an index i into a BS label {n,l} */
  int ntest = 1;
  int ltest = 0;
  int prod = (ntest*(ntest+1))/2;

  while (prod <= i){
   ntest += 1;
   prod = (ntest*(ntest+1))/2; 
   ltest = prod-i-1;
  }

  std::array<int, 2> output;
  output[0] = ntest;
  output[1] = ltest;

  return output;
}

double hypergeometric_prior(double a, double b, double c, double z){
	/* Implementation of eq.(A15) of arXiv:2308.01336v2 */
	const double TOLERANCE = 1.0e-10;
  double term = a*b*z/c;
  double value = 1.0 + term;
  int n = 1;
  while ( abs( term ) > TOLERANCE )
  {
    a++, b++, c++, n++;
    term *= a*b*z/(c*n);
    value += term;
  }
  return value;
}

double* MatrixInversion(int dim, long double *MInput){ 
  /* This function gives the inverse of the asymmetric, real dim*dim matrix MInput using 
  CalcHEP built-in routine rJacobiA, which requires MInput to be of type long double. */
  
  const int total_num_elements = dim*dim;
    // Allocate memory dynamically for MInv
    double *MInvPointer = (double*)malloc(total_num_elements * sizeof(double));
    if (MInvPointer == NULL) {
        fprintf(stderr, "Memory allocation failed for MInvPointer\n");
        return NULL;
    }

    // Allocate memory dynamically for long double arrays
    long double *VMat = (long double*)malloc(total_num_elements * sizeof(long double));
    long double *UMat = (long double*)malloc(total_num_elements * sizeof(long double));
    long double *ValueArray = (long double*)malloc(dim * sizeof(long double));
    if (VMat == NULL || UMat == NULL || ValueArray == NULL) {
        fprintf(stderr, "Memory allocation failed for VMat, UMat, or ValueArray\n");
        free(MInvPointer);
        if (VMat) free(VMat);
        if (UMat) free(UMat);
        if (ValueArray) free(ValueArray);
        return NULL;
    }

    // Initialization
    for (int i = 0; i < dim; i++) {
        ValueArray[i] = 0.;
        for (int j = 0; j < dim; j++) {
            int ind = dim * i + j;
            VMat[ind] = 0.;
            UMat[ind] = 0.;
            MInvPointer[ind] = 0.;
        }
    }

    long double *Vij = VMat;
    long double *Uij = UMat;
    long double *EigVals = ValueArray;

    int sol = rJacobiA(MInput, dim, EigVals, Uij, Vij);
    
    
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            double elem = 0.;
            for (int k = 0; k < dim; k++) {
                elem += (Uij[k * dim + j] * Vij[dim * k + i]) / EigVals[k];
            }
            MInvPointer[dim * i + j] = elem;
        }
    }
    // Free the dynamically allocated memory for intermediate arrays
    free(VMat);
    free(UMat);
    free(ValueArray);
    
  return MInvPointer; //Address of Minv returned 
}

double dofs_particle(char* prtcl_name){ //Returns the degress of freedom of the particle with name "prtcl_name"
  int s1, q1, c1; //Initialise spin and color of the dark sector particle
  qNumbers(prtcl_name, &s1, &q1, &c1); //Assign quantum numbers to dark sector particle)
  return (1. + double(s1))*abs(double(c1)); 
}
