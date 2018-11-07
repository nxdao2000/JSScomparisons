//
//  c code to generate Figure 2 of
//	It is based on source code from Roy Manojit. 
//

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

double expit(double x){ return 1.0/(1.0+exp(-x)); }

double dot_product(int dim, const double *basis, const double *coef)
{
	int j;
	double trans = 0.0;
	for(j=0;j<dim;j++) trans += coef[j]*basis[j];
	return(trans);
}

//////////////////////////////////////////////////////////////
///////////////// measurement models /////////////////////////
//////////////////////////////////////////////////////////////

#define LOG_SIGOBS		(p[parindex[0]])
#define CASE			(x[stateindex[0]])
#define ERCT			(x[stateindex[1]])
#define DATA			(y[obsindex[0]])

void nbinom_rmeasure(double *y, double *x, double *p, int *obsindex, int *stateindex,
			int *parindex, int *covindex, int ncovars, double *covars, double t)
{
	double size, prob, sigOBS;
//	if((ERCT>0.0) || (!(R_FINITE(CASE)))){ DATA = R_NaReal; }	// for debugging
//	else{
		sigOBS = exp(LOG_SIGOBS);		// untransform parameter
		size = 1.0 / (sigOBS * sigOBS);
		prob = size / (CASE + size);
		DATA = rnbinom(size,prob);
//	}
}

void nbinom_dmeasure(double *lik, double *y, double *x, double *p, int give_log,
			int *obsindex, int *stateindex, int *parindex, int *covindex, int ncovars,
			double *covars, double t)
{
	double size, prob, sigOBS;
	static double tol = 1.0e-18;
	if((ERCT>0.0) || (!(R_FINITE(CASE)))){ *lik = tol; } 
	else{
		sigOBS = exp(LOG_SIGOBS);		// untransform parameter
		size = 1.0 / (sigOBS * sigOBS);
		prob = size / (CASE + size);
		*lik = dnbinom(DATA,size,prob+tol,0)+tol;
	}
	if (give_log) *lik = log(*lik);
}

#undef LOG_SIGOBS
#undef CASE
#undef ERCT
#undef DATA


//////////////////////////////////////////////////////////////
//////////////////// process models //////////////////////////
//////////////////////////////////////////////////////////////

/////////////// model SEIH3QS (with three H classes) ////////////////

#define LOG_MUEI		(p[parindex[1]])
#define LOG_MUIH		(p[parindex[2]])
#define LOG_MUHI		(p[parindex[3]])
#define LOG_MUIQ		(p[parindex[4]])
#define LOG_MUIS		(p[parindex[5]])
#define LOG_MUQS		(p[parindex[6]])
#define LOG_TAU			(p[parindex[7]])
#define LOG_SIGPRO		(p[parindex[8]])
#define LOGIT_RHO		(p[parindex[9]])
#define LOGIT_Q			(p[parindex[10]])
#define LOGIT_BR		(p[parindex[11]])
#define DELTA			(p[parindex[12]])
#define B1			(p[parindex[13]])
#define CASE			(x[stateindex[0]])
#define ERCT			(x[stateindex[1]])
#define S			(x[stateindex[2]])
#define E			(x[stateindex[3]])
#define I			(x[stateindex[4]])
#define Q			(x[stateindex[5]])
#define H1			(x[stateindex[6]])
#define H2			(x[stateindex[7]])
#define H3			(x[stateindex[8]])
#define K			(x[stateindex[9]])
#define F			(x[stateindex[10]])
#define W			(x[stateindex[11]])
#define POP			(covar[covindex[0]])
#define DPOPDT			(covar[covindex[1]])
#define RAINFALL		(covar[covindex[2]])
#define SEASON			(covar[covindex[3]])
#define NSPLINE			covdim - 3

void SEIH3QS(double *x, const double *p, const int *stateindex, const int *parindex,
		const int *covindex, int covdim, const double *covar, double t, double dt)
{
	double dW, beta, foi;
	static double muEI, muIH, muHI, muIQ, muIS, muQS, tau, q, br, rho, sigPRO, varPRO,
			delt, rateHI;
	double dBS, dSD, dSE, dEI, dED, dIH1, dIQ, dIS, dID, dH1H2, dH1D, dH2H3, dH2D,
			dH3I, dH3D, dQS, dQD, dK, dF, newcases;

	// untransfomr parameters
	muEI   = exp(LOG_MUEI);
	muIH   = exp(LOG_MUIH);
	muHI   = exp(LOG_MUHI);
	muIQ   = exp(LOG_MUIQ);
	muIS   = exp(LOG_MUIS);
	muQS   = exp(LOG_MUQS);
	tau    = exp(LOG_TAU);
	sigPRO = exp(LOG_SIGPRO);
	rho    = expit(LOGIT_RHO);
	q      = expit(LOGIT_Q);
	br     = expit(LOGIT_BR);

	// compute transition rates
	rateHI = 3*muHI;
	dBS   = DELTA*POP + DPOPDT;
	dSD   = DELTA * S;
	dSE   = F * S;
	dEI   = muEI * E;
	dED   = DELTA * E;
	dIH1  = muIH * I;
	dIQ   = muIQ * I;
	dIS   = muIS * I;
	dID   = DELTA * I;
	dH1H2 = rateHI * H1;
	dH1D  = DELTA * H1;
	dH2H3 = rateHI * H2;
	dH2D  = DELTA * H2;
	dH3I  = rateHI * H3;
	dH3D  = DELTA * H3;
	dQS   = muQS * Q;
	dQD   = DELTA * Q;

	// compute beta
	beta = exp(dot_product(NSPLINE,&SEASON,&B1) + br*RAINFALL);

	// make sure parameters and variables stay bounded
	if(!(R_FINITE(muEI))      || 
		!(R_FINITE(muIH)) ||
		!(R_FINITE(muHI)) ||
		!(R_FINITE(muIQ)) ||
		!(R_FINITE(muIS)) ||
		!(R_FINITE(muQS)) ||
		!(R_FINITE(tau))  ||
		!(R_FINITE(q))    ||      
		!(R_FINITE(br))   ||      
		!(R_FINITE(rho))  ||
		!(R_FINITE(beta)) ||
		!(R_FINITE(S))    ||
		!(R_FINITE(E))    ||
		!(R_FINITE(I))    ||
		!(R_FINITE(Q))    ||
		!(R_FINITE(H1))   ||
		!(R_FINITE(H2))   ||
		!(R_FINITE(H3))   ||
		!(R_FINITE(K))    ||
		!(R_FINITE(F))    ||
		!(R_FINITE(W))    ||
		!(R_FINITE(CASE)))
	return;

	// compute process noise increment
	varPRO  = sigPRO * sigPRO;
	dW = rgamma(dt/varPRO,varPRO);	//gamma noise (mean=dt,variance=varPRO*dt)
	if(!(R_FINITE(dW))) return;

	// compute mosquito force-of-infection
	foi = ((I + q*Q)/POP) * beta * dW/dt;
	delt = 2.0*dt/tau;
	dK = (foi - K)*delt;
	dF = (K - F)*delt;

	// compute equations
	newcases = dEI + dH3I;
	S  += (dBS + dIS + dQS - dSE - dSD)*dt;
	E  += (dSE - dEI - dED)*dt;
	I  += (newcases - dIH1 - dIS - dIQ - dID)*dt;
	H1 += (dIH1 - dH1H2 - dH1D)*dt;
	H2 += (dH1H2 - dH2H3 - dH2D)*dt;
	H3 += (dH2H3 - dH3I - dH3D)*dt;
	Q  += (dIQ - dQS - dQD)*dt;
	K  += dK;
	F  += dF;
	CASE += rho*newcases*dt;
	W += dW;

	// keep track of errors
	if(S<0.0) { ERCT -= S; S=0.0;   }
	if(E<0.0) { ERCT -= E; E=0.0;   }
	if(H1<0.0){ ERCT -= H1; H1=0.0; }
	if(H2<0.0){ ERCT -= H2; H2=0.0; }
	if(H3<0.0){ ERCT -= H3; H3=0.0; }
	if(I<0.0) { ERCT -= I; I=0.0;   }
	if(Q<0.0) { ERCT -= Q; Q=0.0;   }
	if(K<0.0) { ERCT -= K; K=0.0;   }
	if(F<0.0) { ERCT -= F; F=0.0;   }
}

#undef LOG_MUEI
#undef LOG_MUIH
#undef LOG_MUHI
#undef LOG_MUIQ
#undef LOG_MUIS
#undef LOG_MUQS
#undef LOG_TAU
#undef LOG_SIGPRO
#undef LOGIT_RHO
#undef LOGIT_Q
#undef LOGIT_BR
#undef DELTA
#undef B1
#undef CASE
#undef ERCT
#undef S
#undef E
#undef I
#undef Q
#undef K
#undef F
#undef H1
#undef H2
#undef H3
#undef W
#undef POP
#undef DPOPDT
#undef SEASON
#undef RAINFALL
#undef NSPLINE

