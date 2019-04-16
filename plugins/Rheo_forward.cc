/*
 * This module contains the various subroutines
Usage import """
#######################################################################################
*/
//#####################  IMPORT STANDARD MODULES   ######################################   
/*
import sys,os
import numpy as np #for numerical analysis
from scipy import integrate  #integration routine
from subprocess import call #for calling unix commands
from datetime import date  #to give a timestamp to output
import pdb	#for the debugger pdb.set_trace()
from matplotlib.pyplot import * # for the figures
from math import pi,exp,log,sqrt
from ConfigParser import SafeConfigParser
*/
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include <string>
#include <map>
#include <random>

#include "refvalues_and_utilities.h"

using namespace std;



/*
 * ------------------------------------------------------------------
 * 
 * 		CONSTANTS
 * 
 * ------------------------------------------------------------------
 */


// A series of constants relevant for a specific modeling approach.
// The modeling approach name is contained in the name of namespace.
namespace Rheology_Constants
{
	namespace Common
	{
		// Ideal gas constant [J K^-1 mol^-1]
		const double R = 8.3144621;
	}
	//###########################  Jackson & Faul 2010  ###########################
	//###########################   Extended Burgers    ###########################
	namespace JF10_eBurg
	{
		// Shear modulus at TR, PR [Pa]
		const double GUR = 66.5e+9;
		
		// T-derivative of Shear modulus [Pa K^-1]
		const double dGdT = -0.0136e+9;
		
		// P-derivative of Shear modulus
		const double dGdP = 1.8;
		
		// Activation energy ("U" in the paper) [J mol^-1]
		const double E = 360000;
		
		// Activation volume [m^3 mol^-1]
		const double V = 1e-5;
		
		// Reference Temperature [K]
		const double TR = 1173;
		
		// Reference Pressure [Pa]
		const double PR = 0.2e+9;
		
		// Reference Grainsize [m]
		const double gsR = 13.4e-6;
		
		// Anelastic grain size exponent
		const double ma = 1.31;
		
		// Frequency dependence
		const double alpha = 0.274;
		
		// Relaxation strength - Burgers
		const double Delta = 1.04;
		
		// Reference relaxation time - upper bound [s]
		const double tauHR = 1e7;
		
		// Reference relaxation time - lower bound [s]
		const double tauLR = 1e-3;
		
		// Reference relaxation time - Maxwell [s]
		const double tauMR = 3.02e+7;
		
		
		// "extended" dissipation peak parms
		// Peak width
		const double sigma = 4;
		
		// Relaxation strength - Peak
		const double DeltaP = 0.057;
		
		// Reference relaxation time - Peak [s]
		const double tauPR = 3.9811e-4;
		
		// Viscous grain size exponent
		const double mv = 3;
	}
	//###########################  Jackson & Faul 2010  ###########################
	//###########################        Andrade        ###########################
	namespace Andrade
	{
		// Shear modulus at TR, PR [Pa]
		const double GUR = 62.2e+9;
		
		// T-derivative of Shear modulus [Pa K^-1]
		const double dGdT = -0.0136e+9;
		
		// P-derivative of Shear modulus
		const double dGdP = 1.8;
		
		// Activation energy ("U" in the paper) [J mol^-1]
		const double E = 303000;
		
		// Activation volume [m^3 mol^-1]
		const double V = 1e-5;
		
		// Reference Temperature [K]
		const double TR = 1173;
		
		// Reference Pressure [Pa]
		const double PR = 0.2e+9;
		
		// Reference Grainsize [m]
		const double gsR = 3.1e-6;
		
		// grainsize dependence exponent ------------ ?????? who cares ??????
		const double m = 3;
		
		// frequency exponent
		const double n = 0.33;
		
		// Beta star = Beta / J_unrelaxed
		const double Bstar = 0.020;
		
		// Reference relaxation time - Maxwell [s]
		const double tauMR = 1.9953e+5;
	}
	//##############################################################################
	//###########################  McCarthy et al. 2011  ###########################
	namespace M11
	{
		// Shear modulus at TR, PR, Isaak, 1992 [Pa]
		const double GUR = 82e+9;
		
		// T-derivative of Shear modulus, Isaak, 1992 [Pa K^-1]
		const double dGdT = -0.0136e+9;
		
		// P-derivative of Shear modulus, Isaak, 1992
		const double dGdP = 1.8;
		
		// Activation energy ("U" in the paper) [J mol^-1]
		const double E = 505000;
		
		// Activation volume [m^3 mol^-1]
		const double V = 1.2e-5;
		
		// Reference Temperature [K]
		const double TR = 1473;
		
		// Reference Pressure [Pa]
		const double PR = 0.2e+9;
		
		// Reference Grainsize [m]
		const double gsR = 1e-3;
		
		// grainsize dependence exponent
		const double m = 3;
		
		// Reference viscosity [Pa s]
		const double eta0 = 6.6e+19;
		
		// fn polynomial fit parms (equation 26)
		const double a0 = 5.5097e-1;
		const double a1 = 5.4332e-2;
		const double a2 = -2.3615e-3;
		const double a3 = -5.7175e-5;
		const double a4 = 9.9473e-6;
		const double a5 = -3.4761e-7;
		const double a6 = 3.9461e-9;
	}
	//##########################################################################
	//###########################  Takei et al 2014  ###########################
	namespace Tak14
	{
		// Shear modulus at TR, PR [Pa]
		const double GUR = 82e+9;
		
		// T-derivative of Shear modulus [Pa K^-1]
		const double dGdT = -0.0136e+9;
		
		// P-derivative of Shear modulus
		const double dGdP = 1.8;
		
		// Reference Temperature [K]
		const double TR = 1473;
		
		// Reference Pressure (by ref to M11, not in paper??)  Pa
		const double PR = 0.2e+9;
		
		// Reference Grainsize [m]
		const double gsR = 1e-3;
		
		// Activation energy ("H" in the paper) [J mol^-1]
		const double E = 505000;
		
		// Activation volume [m^3 mol^-1]
		const double V = 1.2e-5;
		
		// Reference Viscosity [Pa s]
		const double eta0 = 6.6e+19;
		
		// grainsize dependence exponent
		const double m = 3;
		
		// Peak standard deviation (value for melt-free Ol)
		const double sigmap = 4;
		
		// Peak pre-exponent (value for melt-free Ol)
		const double Ap = 0.007;
	}
	//###################################################################################
	//###########################  Priestley & McKenzie 2013  ###########################
	namespace P_M13
	{
		// Shear modulus at TR, PR [Pa]
		const double GUR = 72.66e+9;
		
		// T-derivative of Shear modulus [Pa K^-1]
		const double dGdT = -0.00871e+9;
		
		// P-derivative of Shear modulus
		const double dGdP = 2.04;
		
		// Activation energy ("U" in the paper) [J mol^-1]
		const double E = 402900;
		
		// Activation volume - m^3 mol^-1]
		const double V = 7.81e-6;
		
		// Reference Temperature [K]
		const double TR = 1473;
		
		// Reference Pressure [Pa]
		const double PR = 1.5e+9;
		
		// Reference Viscosity [Pa s]
		const double eta0 = 2.3988e+22;
		
		// Reference Grainsize [m]
		const double gsR = 1e-3;
		
		// grainsize dependence exponent
		const double m = 3;
		
		// fn polynomial fit parms (equation 26)
		const double a0 = 5.5097e-1;
		const double a1 = 5.4332e-2;
		const double a2 = -2.3615e-3;
		const double a3 = -5.7175e-5;
		const double a4 = 9.9473e-6;
		const double a5 = -3.4761e-7;
		const double a6 = 3.9461e-9;
	}
	//##############################################################################
	//###########################  LOWER MANTLE GUESSES  ###########################
	namespace LOWM
	{
		// Activation volume - m^3 mol^-1
		const double V = 1e-6;
	}
}


/*
 * ------------------------------------------------------------------
 * 
 *		MODULES
 * 
 *-------------------------------------------------------------------
 */


std::pair<double, double> PTd_VsQ( double freq, double P, double T, double gs, double rho, std::string model, double highQapprox, double Vs0 = std::numeric_limits<double>::quiet_NaN() )
{
	/*
	 * This subroutine calculates the values of Vs and Q at a given 
	 * PTd condition at a frequency f and density rho. If highQapprox is 1
	 * approximate values of Vs and Q can be calculated from 
	 * e.g. McCarthy et al. (2011) eqn. 24, else use their more explicit form eqn. B6
	 */
	
	//---------[ NB this uses V_activation = sloping from 11 - 3 (e-6) in the LM ]---------
	
	
	double J1byJu, J2byJu, Ju, Vs, G;
		
	// Value of the real and complex part of compliance			
	if ( std::isnan(Vs0) )
	{
		J1byJu = getJ1byJu(freq, P, T, gs, model);
		J2byJu = getJ2byJu(freq, P, T, gs, model);
		Ju = getJu(T, P, model);
	}
	else
	{
		J1byJu = getJ1byJu(freq, P, T, gs, model, Vs0, rho);
		J2byJu = getJ2byJu(freq, P, T, gs, model, Vs0, rho);
		Ju = getJu(T, P, model, Vs0, rho);
	}
	
	
	double J1 = J1byJu * Ju;
	double J2 = J2byJu * Ju;
	
	
	// EVERYONE AGREES
	double qs = J2 / J1;
	
	
	// VS DISAGREEMENT	
	if (model == "M11") Vs = 1.0 / ( sqrt(rho * J1) );
	else if (model == "P_M13") Vs = 1.0/( sqrt(rho * J1) );
	else if (model == "JF10_eBurg")
	{
		G = pow((J1*J1 + J2*J2), -0.5);
		Vs = sqrt(G/rho);
	}
	else if (model == "Tak14")
	{
		// pdb.set_trace()
		G = pow(J1*J1 + J2*J2, -0.5);
		Vs = sqrt(G/rho);
	}
	
	
	// **** When Qs^-1 << 1 (J2/J1 << 1), highQapprox can be 1. If Qs^-1 is not small, Qs^-1 does not equal Qmu^-1 ****
	if (highQapprox == 0)
	{
		double factor = ( 1.0 + sqrt(1.0 + pow((J2/J1), 2.0) ) )/2.0;
		qs = qs/factor;
		Vs = Vs/sqrt(factor);
	}
	
	
	double Q = 1.0 / qs;
	
	
	// if (Q > 1e9) pdb.set_trace();
	
	return std::make_pair(Vs,Q);
}



double getJ1byJu(double freq, double P, double T, double gs, std::string model, double Vs0 = std::numeric_limits<double>::quiet_NaN(), double rho = std::numeric_limits<double>::quiet_NaN() )
{
	// This subroutine gets the value of J1/Ju based on expressions in a paper
	
	
	double J1byJu, Ju, omega;
	double fn, val;
	double a0, a1, a2, a3, a4, a5, a6;
	
	
	// If model is McCarthy 2011
	if (model == "M11")
	{
		if ( std::isnan(Vs0) ) fn = getfn(freq, P, T, gs, model);
		else fn = getfn(freq, P, T, gs, model, Vs0, rho);
		
		
		// Getting J1 values from McCarthy et al., 2011
		a0 = getconstant(model, "a0");
		a1 = getconstant(model, "a1");
		a2 = getconstant(model, "a2");
		a3 = getconstant(model, "a3");
		a4 = getconstant(model, "a4");
		a5 = getconstant(model, "a5");
		a6 = getconstant(model, "a6");

		if (fn <= 1.0e+13)
		{
			val = 0.0;
			val = val + a0*(pow(log(fn), 0));
			val = val + a1*(pow(log(fn), 1));
			val = val + a2*(pow(log(fn), 2));
			val = val + a3*(pow(log(fn), 3));
			val = val + a4*(pow(log(fn), 4));
			val = val + a5*(pow(log(fn), 5));
			val = val + a6*(pow(log(fn), 6));
			J1byJu = 1.0 / val;
		}
		else J1byJu = 1.0;
	}
//#	If model is Priestley and McKenzie 2013		
	else if (model == "P_M13")
	{
		fn = getfn(freq, P, T, gs, model, Vs0, rho);
		
		// Getting J1 values from McCarthy et al., 2011
		a0 = getconstant(model, "a0");
		a1 = getconstant(model, "a1");
		a2 = getconstant(model, "a2");
		a3 = getconstant(model, "a3");
		a4 = getconstant(model, "a4");
		a5 = getconstant(model, "a5");
		a6 = getconstant(model, "a6");
		
		if (fn <= 1.0e+13)
		{
			val = 0.0;
			val = val + a0*(pow(log(fn), 0));
			val = val + a1*(pow(log(fn), 1));
			val = val + a2*(pow(log(fn), 2));
			val = val + a3*(pow(log(fn), 3));
			val = val + a4*(pow(log(fn), 4));
			val = val + a5*(pow(log(fn), 5));
			val = val + a6*(pow(log(fn), 6));
			J1byJu = 1.0 / val;
		}
		else J1byJu = 1.0;
	}
//#	If model is extended Burgers from Jackson & Faul 2010
	else if (model == "JF10_eBurg")
	{
		double alpha = getconstant(model,"alpha");
		double sigma = getconstant(model,"sigma");
		double Delta = getconstant(model,"Delta");
		double tauHR = getconstant(model,"tauHR");
		double tauLR = getconstant(model,"tauLR");
		double tauPR = getconstant(model,"tauPR");
		double DeltaP = getconstant(model,"DeltaP");
		double ma = getconstant(model,"ma");
		
		omega = 2.0 * UniversalConst::PI * freq;
		double tauH = gettau(tauHR, ma, P, T, gs, model);
		double tauL = gettau(tauLR, ma, P, T, gs, model);
		double tauP = gettau(tauPR, ma, P, T, gs, model);
		
		double j1b = (alpha*Delta)/(pow(tauH,alpha) - pow(tauL,alpha));
		double j1p = DeltaP/(sigma*sqrt(2*UniversalConst::PI));
		
		// !!!!! REPLACE THESE !!!!!
		double i1b = intgrate_j1b(tauL,tauH,alpha,omega,N=1000);
		double i1p = integrate.quad(intgrnd_j1p, 0, np.inf, args=(omega,tauP,sigma))[0];
		
		
		J1byJu = 1.0 + (j1b*i1b) + (j1p*i1p);
	}
//#	If model is Takei 2014
	else if (model == "Tak14")
	{
		double m = getconstant(model, "m"); // grainsize exponent
		double eta0 = getconstant(model, "eta0"); // ref viscosity
		double Ap = getconstant(model, "Ap"); // pre-exponent for peak (normalisation factor)
		double sigmap = getconstant(model, "sigmap"); // controls peak width - see Takei 2014 eqn 17
		
		// Value of the real and complex part of compliance			
		if ( std::isnan(Vs0) ) Ju = getJu(T, P, model);
		else Ju = getJu(T, P, model, Vs0, rho);
		
		omega = 2.0 * UniversalConst::PI * freq;
		double tauM = gettau(Ju*eta0, m, P, T, gs, model);
		
		double i1 = integrate.quad(intgrnd_j1tak, 0.0, np.inf, args=(omega, tauM, Ap, sigmap))[0];

		J1byJu = 1.0 + i1;
	}
	
	
	return J1byJu;
}



double getJ2byJu(double freq, double P, double T, double gs, std::string model, double Vs0 = std::numeric_limits<double>::quiet_NaN(), double rho = std::numeric_limits<double>::quiet_NaN() )
{
	// This subroutine gets the value of J2 based on expressions in a paper
	
	double fn, taun, Xn, J2byJu, omega, tauM;
	
	if (model == "M11")  // *** If model is McCarthy 2011
	{
		fn = getfn(freq, P, T, gs, model, Vs0, rho);
		taun = 1.0/(2.0 * UniversalConst::PI * fn);
		
		// Getting J2 values from McCarthy et al., 2011
		Xn = mccarthyXn(taun);
		J2byJu = (UniversalConst::PI / 2.0*Xn) + taun / (2.0 * UniversalConst::PI);
	}
	else if (model == "P_M13")  // *** If model is Priestley and McKenzie 2013
	{
		fn = getfn(freq, P, T, gs, model, Vs0, rho);
		taun = 1.0/(2.0 * UniversalConst::PI * fn);
		
		//#Getting J2 values from McCarthy et al., 2011
		Xn = mccarthyXn(taun);
		J2byJu = (UniversalConst::PI/2.0*Xn)+taun/(2.0*UniversalConst::PI);
	}
	else if (model == "JF10_eBurg")  // *** If model is extended Burgers from Jackson & Faul 2010
	{
		double alpha = getconstant(model, "alpha");
		double sigma = getconstant(model, "sigma");
		double Delta = getconstant(model, "Delta");
		double tauHR = getconstant(model, "tauHR");
		double tauLR = getconstant(model, "tauLR");
		double tauPR = getconstant(model, "tauPR");
		double tauMR = getconstant(model, "tauMR");
		double DeltaP = getconstant(model, "DeltaP");
		double ma = getconstant(model, "ma");
		double mv = getconstant(model, "mv");
		
		omega = 2.0 * UniversalConst::PI * freq;
		
		double tauH = gettau(tauHR, ma, P, T, gs, model);
		double tauL = gettau(tauLR, ma, P, T, gs, model);
		double tauP = gettau(tauPR, ma, P, T, gs, model);
		tauM = gettau(tauMR, mv, P, T, gs, model);
		
		double j2b = omega*(alpha*Delta)/(pow(tauH,alpha) - pow(tauL,alpha));
		double j2p = omega*DeltaP/(sigma*sqrt(2*UniversalConst::PI));
		
		// !!!!! NEED TO REPLACE ITEMS HERE !!!!!
		std::vector<double> tau_arr = MyUtilities::linspace(tauL,tauH,100); //tau_arr=np.linspace(tauL,tauH,100);
		double i2b = intgrate_j2b(tauL, tauH, alpha, omega, int N=1000);
		double i2p = integrate.quad(intgrnd_j2p, 0, np.inf, args=(omega,tauP,sigma))[0];
		
		J2byJu = (j2b*i2b)+(j2p*i2p)+(1.0/(omega*tauM));
	}
	else if (model == "Tak14")  // *** If model is Takei 2014
	{
		double m = getconstant(model,"m"); // grainsize exponent
		double eta0 = getconstant(model,"eta0"); // ref viscosity
		double Ap = getconstant(model,"Ap"); // pre-exponent for peak (normalisation factor)
		double sigmap = getconstant(model,"sigmap"); // controls peak width - see Takei 2014 eqn 17
		
		
		// Value of the real and complex part of compliance			
		double Ju;
		if ( std::isnan(Vs0) ) Ju = getJu(T, P, model);
		else Ju = 1.0 / (rho * pow(Vs0, 2.0));
		
		
		omega = 2.0 * UniversalConst::PI * freq;
		tauM = gettau(Ju*eta0, m, P, T, gs, model);
		
		// !!!!! NEED TO REPLACE ITEMS HERE !!!!!
		double i2 = integrate.quad(intgrnd_j2tak, 0, np.inf, args=(omega,tauM,Ap,sigmap))[0];
		
		
		J2byJu = i2 + (1.0 / (omega*tauM) );
	}
	
	
	return J2byJu;
}



double mccarthyXn(double taun)
{
	// This subroutine gets value of relaxation spectrum at a value of 
	// normalized time scale (taun) from McCarthy et al. (2011) eqn. 25
	
	double Xn;
	
	if (taun >= 1.0e-11) Xn = 0.32 * (pow(taun, (0.39-0.28/(1.0+2.6 * pow(taun, 0.1) ) ) ) );
	else Xn = 1853.0 * sqrt(taun);
	
	
	return Xn;
}



double getfn(double freq, double P, double T, double gs, std::string model, double Vs0 = std::numeric_limits<double>::quiet_NaN(), double rho = std::numeric_limits<double>::quiet_NaN() )
{
	// Get the normalized frequency from McCarthy et al. (2011) eqn. 19 and 22
	
	double eta0 = getconstant(model,"eta0");
	double m = getconstant(model,"m");
	
	double Ju;
	
	if ( std::isnan(Vs0) ) Ju = getJu(T,P,model);
	else Ju = getJu(T,P,model,Vs0,rho);
	
	
	double tau0 = Ju * eta0;
	double taur = gettau(tau0, m, P, T, gs, model);
	
	if (model == "P_M13") taur = gettauM_P_M13(T,P,Ju);
	
	double fn = freq * taur;
	
	return fn;
}



double geteta(double P, double T, double gs, std::string model)
{
	// Get the viscosity from a combination of power law and Arrhenius eqn. like 
	// in McCarthy eqn. 8
	
	double TR = getconstant(model,"TR");
	double gsR = getconstant(model,"gsR");
	double eta0 = getconstant(model,"eta0");
	
	if (model == "M11" || model == "P_M13") m=getconstant(model,"m");
	else if (model == "eBurgers") m=getconstant(model,"ma");
	
	double E = getconstant(model,"E");
	double R = getconstant("Common","R");
	
	double eta = eta0*pow((gs/gsR),m)*exp((E/R)*(TR-T)/(TR*T));
	
	return eta;
}



double getJu(double T, double P, std::string model, double Vs0 = std::numeric_limits<double>::quiet_NaN(), double rho = std::numeric_limits<double>::quiet_NaN() )
{
	// This gets the unrelaxed modulus at T,P conditions based on constants from study
	
	double GUR, dGdT, dGdP, TR, PR;
	double Ju;
	
	if ( std::isnan(Vs0) )
	{
		GUR = getconstant(model,"GUR");
		dGdT = getconstant(model,"dGdT"); 
		dGdP = getconstant(model,"dGdP");
		TR = getconstant(model,"TR");
		PR = getconstant(model,"PR");
		Ju = 1.0*(GUR + dGdT*(T-TR) + dGdP*(P-PR));
	}
	else Ju = 1.0/(rho*pow(Vs0, 2.0));
	
	return Ju;
}



double gettau(tau0, m, double P, double T, double gs, std::string model)
{
	// calculate the relaxation time at P,T,gs, given model and reference(0) relaxation 
	// time and grainsize exponent
	double PR = getconstant(model,"PR");
	double TR = getconstant(model,"TR");
	double gsR = getconstant(model,"gsR");
	double E = getconstant(model,"E");
	double V = getconstant(model,"V");
	double R = getconstant("Common","R");
	
	if (P > 24.3e9)
	{
		V = 5e-6;	//# constant
//# 		V = -7.2727e-17*P + 1.1770e-05		# 10 to 2
//# 		V = -5.4545e-17*P + 1.13250e-05		# 10 to 4
//# 		V = -7.2727e-17*P + 1.3767e-05		# 12 to 4
//# 		Vact_lowM(P,Vtop,Vbot,Ptop,Pfold) where Pfold of 11e9 gives 1/2 LM die-out
//# 		V = Vact_lowM(P,15e-6,4e-6,24.3e9,11e9)
	}
	
	double tau = tau0 * pow((gs/gsR), m) * exp( (E/R)*(TR-T)/(TR*T) + (V/R)*(P*TR-PR*T)/(T*TR) );
	
	return tau;
}



double gettauM_P_M13(double T, double P, double Ju)  // !!!!! getconstant needs replacement. See next function !!!!!
{
	// This gets the Maxwell relaxation time at T,P conditions based on constants
	// from study model='P_M13'
	double TR = getconstant(model,"TR");
	double PR = getconstant(model,"PR");
	double E = getconstant(model,"E");
	double V = getconstant(model,"V");
	double R = getconstant("Common","R");
	double eta0 = getconstant(model,"eta0");
	
	
	if (P > 24.3e9)
	{
		V = 5e-6;	//# constant
//# 		V = -7.2727e-17*P + 1.1770e-05		# 10 to 2
//# 		V = -5.4545e-17*P + 1.13250e-05		# 10 to 4
//# 		V = -7.2727e-17*P + 1.3767e-05		# 12 to 4
//# 		Vact_lowM(P,Vtop,Vbot,Ptop,Pfold) where Pfold of 11e9 gives 1/2 LM die-out
//# 		V = Vact_lowM(P,15e-6,4e-6,24.3e9,11e9)
	}
	
	
	double astar = exp((E - PR*V)/(R*TR))/exp((E-P*V)/(R*T));
	double eta = eta0/astar;
	double tauM = Ju*eta;
	
	return tauM;
}


/*
double getconstant(std::string model, std::string par)  // !!!!! THIS FUNCTION NEEDS TO BE REPLACED. NO SafeConfigParser in c++ !!!!!
{
	// This subroutine gets the value of a constant relevant for a modeling approach
	parser = SafeConfigParser();
	parser.read("constants_complete.ini");
	val=float(parser.get(model,par));
	
	return val;
}
*/


/*
 * ------------------------------------------------------------------
 * 
 *		Integrations for JF10
 * 
 *-------------------------------------------------------------------
 */

/**
 * Function to calculate the function from J1B at a
 * specific coordinate xx.
 * 
 * Needs to be integrated over xx
 */
template <class T>
struct Function_J1B
{
	double alpha, omega;
	
	Function_J1B (double A, double w) : alpha(A), omega(w) {}
	double operator()(double xx)
	{
		return ( pow(xx, alpha-1.0) / (1.0 + pow(omega*xx, 2.0) ) );
	}
};
double intgrate_j1b(double limL, double limH, double alpha, double omega, int N=1000) //*
{
	std::vector<double xx = MyUtilities::logspace(log10(limL), log10(limH), N, 10.0);
	std::vector<double> yy(xx.size(),0.0);
	//for (unsigned int ii=0; ii<xx.size(); ++ii)
	//{
	//	yy[ii] = pow(xx[ii], alpha-1.0) / (1.0 + pow(omega*xx[ii], 2.0) );
	//}
	
	
	Function_J1B<double> J1Bfunc(alpha, omega);
	
	//yy[0] = 0.0; // !!!!! The vector is already initialized to zero, maybe this can go !!!!!
	
	
	// !!! SHIFT TO START FROM 1 TO XX.SIZE() AND FIRST ELEMENT IS 0 !!!
	for (unsigned int inc=1; inc<xx.size(); ++inc)
	{
		unsigned int inc_prev = inc - 1;
		double inc1 = xx[inc_prev]; //radius[irad];
		double inc2 = xx[inc]; //radius[irad+1];
		
		
		// Calculate the cumulative integral
		double current_increment = IntegrationRoutines::qromb(J1Bfunc,inc1,inc2, 1.0e-10);
		if (inc == 0)
		{
			yy[inc] = current_increment;
		}
		else
		{
			yy[inc] = yy[inc-1] + current_increment;
		}
	}
	
	// lines are from the python script. These can be deleted when everything else tested
	//I=integrate.cumtrapz(yy,xx)
	//I = I[len(I)-1]
	//return I;
	
	return yy[xx.size()-1]; // Only return the last element
}



/**
 * Function to calculate the function from J2B at a
 * specific coordinate xx.
 * 
 * Needs to be integrated over xx
 */
template <class T>
struct Function_J2B
{
	double alpha, omega;
	
	Function_J2B (double A, double w) : alpha(A), omega(w) {}
	double operator()(double xx)
	{
		return ( pow(xx, alpha) / (1.0 + pow(omega*xx, 2.0) ) );
	}
};
double intgrate_j2b(double limL, double limH, double alpha, double omega, int N=1000) //*
{
	std::vector<double> xx = MyUtilities::logspace(log10(limL), log10(limH), N, 10.0);
	std::vector<double> yy(xx.size(),0.0);
	/*
	for (unsigned int ii=0; ii<xx.size(); ++ii)
	{
		yy[ii] = pow(xx[ii], alpha) / (1.0 + pow(omega*xx[ii], 2.0) );
	}
	*/
	
	Function_J2B<double> J2Bfunc(alpha, omega);
	
	
	// !!! SHIFT TO START FROM 1 TO XX.SIZE() AND FIRST ELEMENT IS 0 !!!
	for (unsigned int inc=1; inc<xx.size(); ++inc)
	{
		unsigned int inc_prev = inc - 1;
		double inc1 = xx[inc_prev];
		double inc2 = xx[inc];
		
		
		// Calculate the cumulative integral
		double current_increment = IntegrationRoutines::qromb(J2Bfunc,inc1,inc2, 1.0e-10);
		if (inc == 0)
		{
			yy[inc] = current_increment;
		}
		else
		{
			yy[inc] = yy[inc-1] + current_increment;
		}
	}
	/*
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I;
	*/
	return yy[xx.size()-1]; // Only return the last element
}



/**
 * Function to calculate the function from J1P at a
 * specific coordinate xx.
 * 
 * Needs to be integrated over xx
 */
template <class T>
struct Function_J1P
{
	double omega, tauP, sigma;
	
	Function_J1P (double w, double Tau, double Sig) : omega(w), tauP(Tau), sigma(Sig) {}
	double operator()(double xx)
	{
		return ( (1.0/xx) * (1.0/(1.0 + pow(omega*xx,2.0))) * exp(-pow(log(xx/tauP), 2.0) / (2.0*pow(sigma, 2.0))) );
	}
};
double intgrate_j1p(double limL, double limH, double omega, double tauP, double sigma, int N=1000) //*
{
	std::vector<double> xx = MyUtilities::logspace(log10(limL), log10(limH), N, 10.0);
	std::vector<double> yy(xx.size(),0.0);
	/*
	for (unsigned int ii=0; ii<xx.size(); ++ii)
	{
		yy[ii] = (1.0/xx[ii]) * (1.0/(1.0 + pow(omega*xx[ii],2.0))) * exp(-pow(log(xx[ii]/tauP), 2.0) / (2.0*pow(sigma, 2.0)));
	}
	*/
	
	Function_J1P<double> J1Pfunc(omega, tauP, sigma);
	
	
	// !!! SHIFT TO START FROM 1 TO XX.SIZE() AND FIRST ELEMENT IS 0 !!!
	for (unsigned int inc=1; inc<xx.size(); ++inc)
	{
		unsigned int inc_prev = inc - 1;
		double inc1 = xx[inc_prev];
		double inc2 = xx[inc];
		
		
		// Calculate the cumulative integral
		double current_increment = IntegrationRoutines::qromb(J1Pfunc,inc1,inc2, 1.0e-10);
		if (inc == 0)
		{
			yy[inc] = current_increment;
		}
		else
		{
			yy[inc] = yy[inc-1] + current_increment;
		}
	}
	
	/*
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I;
	*/
	return yy[xx.size()-1]; // Only return the last element
}



 /**
 * Function to calculate the function from J2P at a
 * specific coordinate xx.
 * 
 * Needs to be integrated over xx
 */
/*
template <class T>
struct Function_J2P
{
	double omega, tauP, sigma;
	
	Function_J2P (double w, double Tau, double Sig) : omega(w), tauP(Tau), sigma(Sig) {}
	double operator()(double xx)
	{
		return ( (1.0/(1.0 + pow(omega*xx,2.0)))*exp(-pow(log(xx/tauP),2.0)/(2.0*pow(sigma,2.0))) );
	}
}; 
double intgrate_j2p(double limL, double limH, double omega, double tauP, double sigma, int N=1000) //*
{
	std::vector<double> xx = MyUtilities::logspace(log10(limL), log10(limH), N, 10.0);
	std::vector<double> yy(xx.size(),0.0);
	
	//for (unsigned int ii=0; ii<xx.size(); ++ii)
	//{
	//	yy[ii] = (1.0/(1.0 + pow(omega*xx[ii],2.0)))*exp(-pow(log(xx[ii]/tauP),2.0)/(2.0*pow(sigma,2.0)));
	//}
	
	
	
	Function_J2P<double> J2Pfunc(omega, tauP, sigma);
	
	
	// !!! SHIFT TO START FROM 1 TO XX.SIZE() AND FIRST ELEMENT IS 0 !!!
	for (unsigned int inc=1; inc<xx.size(); ++inc)
	{
		unsigned int inc_prev = inc - 1;
		double inc1 = xx[inc_prev];
		double inc2 = xx[inc];
		
		
		// Calculate the cumulative integral
		double current_increment = IntegrationRoutines::qromb(J2Pfunc,inc1,inc2, 1.0e-10);
		if (inc == 0)
		{
			yy[inc] = current_increment;
		}
		else
		{
			yy[inc] = yy[inc-1] + current_increment;
		}
	}
	
	//I=integrate.cumtrapz(yy,xx)
	//I = I[len(I)-1]
	//return I;
	
	return yy[xx.size()-1]; // Only return the last element
}
*/



double intgrnd_j1b(double tau, double alpha, double omega) //*
{
	return ( ( pow(tau,alpha-1.0) )/(1.0 + pow( ( (omega*tau), 2.0) ) ) );
}



double intgrnd_j1p(double tau, double omega, double tauP, double sigma) //*
{
	return ( (1.0/tau) * (1.0 / (1.0 + pow(omega*tau, 2.0) ) ) * exp(-pow(log(tau/tauP), 2.0) / (2.0*pow(sigma, 2.0) ) ) );
}



double intgrnd_j2b(double tau, double alpha, double omega) //*
{
	return ( pow(tau, alpha) / (1.0 + pow(omega*tau, 2.0) ) );
}



double intgrnd_j2p(double tau, double omega, double tauP, double sigma) //*
{
	return ( (1.0/(1.0 + pow(omega*tau, 2.0) ) ) * exp(-pow(log(tau/tauP), 2.0) / (2.0*pow(sigma, 2.0) ) ) );
}



/*
 * ------------------------------------------------------------------
 * 
 *		Integrations for Takei 2014
 * 
 *-------------------------------------------------------------------
 */

double intgrnd_j1tak(double tau, double omega, double tauM, double Ap, double sigmap) //*
{
	// Integration part of J1byJu term (eqn 8 in Takei et al 2014)
	return ( (1.0 / tau)*(1.0 / (1.0 + pow(omega*tau, 2.0))) * (0.444*pow(tau/tauM, 0.38) + Ap*exp(-0.5*pow(((log(tau/tauM) + 8.1)/sigmap), 2.0))) );
}



double intgrnd_j2tak(double tau, double omega, double tauM, double Ap, double sigmap) //*
{
	// Integration part of J2byJu term (eqn 8 in Takei et al 2014)
	return ( (omega/(1.0 + pow(omega*tau, 2.0) ) ) * (0.444*pow(tau/tauM, 0.38) + Ap*exp(-0.5*pow(((log(tau/tauM) + 8.1) / sigmap), 2.0) ) ) );
}



double Vact_lowM(double P, double Vtop, double Vbot, double Ptop, double Pfold) //*
{
	return ( Vbot + (Vtop-Vbot) * exp( (Ptop-P) / Pfold ) );
}




/*
 ************************************************************
 * 
 * 
 * A function to fetch the correct constant value for the
 * specified Model and Parameter.
 * 
 * !!!!!
 * 			The constants are simple copy and pasted here
 * 			from the constants_complete.ini file.
 * 			These values should not be hardcoded, and
 * 			should be changed later to read from file.
 * !!!!!
 * 
 ************************************************************
 */

double getconstant(std::string model, std::string par)  // !!!!! THIS FUNCTION NEEDS TO BE REPLACED. NO SafeConfigParser in c++ !!!!!
{
	// This subroutine gets the value of a constant relevant for a modeling approach
	
	//val=float(parser.get(model,par));
	double val;
	
	if (model == "Common")
	{
		// Ideal gas constant							J K^-1 mol^-1
		if (par == "R") val = 8.3144621;
	}
	//###########################  Jackson & Faul 2010  ###########################
	//###########################   Extended Burgers    ###########################
	else if (model == "JF10_eBurg")
	{
		//# Shear modulus at TR, PR 						Pa
		if (par == "GUR") val = 66.5e+9;

		//# T-derivative of Shear modulus					Pa K^-1
		else if (par == "dGdT") val = -0.0136e+9;

		//# P-derivative of Shear modulus					
		else if (par == "dGdP") val = 1.8;

		//# Activation energy	("U" in the paper)			J mol^-1
		else if (par == "E") val = 360000;

		//# Activation volume								m^3 mol^-1
		else if (par == "V") val = 1e-5;

		//# Reference Temperature							K
		else if (par == "TR") val = 1173;

		//# Reference Pressure							Pa
		else if (par == "PR") val = 0.2e+9;

		//# Reference Grainsize							m
		else if (par == "gsR") val = 13.4e-6;

		//# Anelastic grain size exponent	
		else if (par == "ma") val = 1.31;

		//# Frequency dependence	
		else if (par == "alpha") val = 0.274;

		//# Relaxation strength - Burgers					
		else if (par == "Delta") val = 1.04;

		//# Reference relaxation time - upper bound		s
		else if (par == "tauHR") val = 1e7;

		//# Reference relaxation time - lower bound		s
		else if (par == "tauLR") val = 1e-3;

		//# Reference relaxation time - Maxwell			s
		else if (par == "tauMR") val = 3.02e+7;


		//# "extended" dissipation peak parms
		//# Peak width
		else if (par == "sigma") val = 4;

		//# Relaxation strength - Peak 					
		else if (par == "DeltaP") val = 0.057;

		//# Reference relaxation time - Peak			s
		else if (par == "tauPR") val = 3.9811e-4;

		//# Viscous grain size exponent
		else if (par == "mv") val = 3;
	}
	//###########################  Jackson & Faul 2010  ###########################
	//###########################        Andrade        ###########################
	else if (model == "Andrade")
	{
		//# Shear modulus at TR, PR 						Pa
		if (par == "GUR") val = 62.2e+9;

		//# T-derivative of Shear modulus					Pa K^-1
		else if (par == "dGdT") val = -0.0136e+9;

		//# P-derivative of Shear modulus					
		else if (par == "dGdP") val = 1.8;

		//# Activation energy	("U" in the paper)			J mol^-1
		else if (par == "E") val = 303000;

		//# Activation volume								m^3 mol^-1
		else if (par == "V") val = 1e-5;

		//# Reference Temperature							K
		else if (par == "TR") val = 1173;

		//# Reference Pressure							Pa
		else if (par == "PR") val = 0.2e+9;

		//# Reference Grainsize							m
		else if (par == "gsR") val = 3.1e-6;

		//# grainsize dependence exponent 		??????	 who cares 	??????
		else if (par == "m") val = 3;

		//# frequency exponent
		else if (par == "n") val = 0.33;

		//# Beta star = Beta / J_unrelaxed	
		else if (par == "Bstar") val = 0.020;

		//# Reference relaxation time - Maxwell			s
		else if (par == "tauMR") val = 1.9953e+5;
	}

	//###########################  McCarthy et al. 2011  ###########################
	else if (model == "M11")
	{
		//# Shear modulus at TR, PR, Isaak, 1992			Pa
		if (par == "GUR") val = 82e+9;

		//# T-derivative of Shear modulus,  Isaak, 1992	Pa K^-1
		else if (par == "dGdT") val = -0.0136e+9;

		//# P-derivative of Shear modulus,  Isaak, 1992	
		else if (par == "dGdP") val = 1.8;

		//# Activation energy	("U" in the paper)			J mol^-1
		else if (par == "E") val = 505000;

		//# Activation volume								m^3 mol^-1
		else if (par == "V") val = 1.2e-5;

		//# Reference Temperature							K
		else if (par == "TR") val = 1473;

		//# Reference Pressure							Pa 
		else if (par == "PR") val = 0.2e+9;

		//# Reference Grainsize							m
		else if (par == "gsR") val = 1e-3;

		//# grainsize dependence exponent
		else if (par == "m") val = 3;

		//# Reference viscosity							Pa s
		else if (par == "eta0") val = 6.6e+19;

		//# fn polynomial fit parms (equation 26)
		else if (par == "a0") val = 5.5097e-1;
		else if (par == "a1") val = 5.4332e-2;
		else if (par == "a2") val = -2.3615e-3;
		else if (par == "a3") val = -5.7175e-5;
		else if (par == "a4") val = 9.9473e-6;
		else if (par == "a5") val = -3.4761e-7;
		else if (par == "a6") val = 3.9461e-9;
	}

	//###########################  Takei et al 2014  ###########################
	else if (model == "Tak14")
	{
		//# Shear modulus at TR, PR 						Pa
		if (par == "GUR") val = 82e+9;

		//# T-derivative of Shear modulus					Pa K^-1
		else if (par == "dGdT") val = -0.0136e+9;

		//# P-derivative of Shear modulus					
		else if (par == "dGdP") val = 1.8;

		//# Reference Temperature							K
		else if (par == "TR") val = 1473;

		//# Reference Pressure (by ref to M11, not in paper??)  Pa
		else if (par == "PR") val = 0.2e+9;

		//# Reference Grainsize							m
		else if (par == "gsR") val = 1e-3;

		//# Activation energy	("H" in the paper)			J mol^-1
		else if (par == "E") val = 505000;

		//# Activation volume								m^3 mol^-1
		else if (par == "V") val = 1.2e-5;

		//# Reference Viscosity							Pa s
		else if (par == "eta0") val = 6.6e+19;

		//# grainsize dependence exponent 		
		else if (par == "m") val = 3;

		//# Peak standard deviation (value for melt-free Ol)
		else if (par == "sigmap") val = 4;

		//# Peak pre-exponent	(value for melt-free Ol)
		else if (par == "Ap") val = 0.007;
	}


	//###########################  Priestley & McKenzie 2013  ###########################
	else if (model == "P_M13")
	{
		//# Shear modulus at TR, PR, 						Pa
		if (par == "GUR") val = 72.66e+9;

		//# T-derivative of Shear modulus,  				Pa K^-1
		else if (par == "dGdT") val = -0.00871e+9;

		//# P-derivative of Shear modulus,  				
		else if (par == "dGdP") val = 2.04;

		//# Activation energy	("U" in the paper)			J mol^-1
		else if (par == "E") val = 402900;

		//# Activation volume								m^3 mol^-1
		else if (par == "V") val = 7.81e-6;

		//# Reference Temperature							K
		else if (par == "TR") val = 1473;

		//# Reference Pressure							Pa 
		else if (par == "PR") val = 1.5e+9;

		//# Reference Viscosity							Pa s
		else if (par == "eta0") val = 2.3988e+22;

		//# Reference Grainsize							m
		else if (par == "gsR") val = 1e-3;

		//# grainsize dependence exponent
		else if (par == "m") val = 3;

		//# fn polynomial fit parms (equation 26)
		else if (par == "a0") val = 5.5097e-1;
		else if (par == "a1") val = 5.4332e-2;
		else if (par == "a2") val = -2.3615e-3;
		else if (par == "a3") val = -5.7175e-5;
		else if (par == "a4") val = 9.9473e-6;
		else if (par == "a5") val = -3.4761e-7;
		else if (par == "a6") val = 3.9461e-9;
	}

	//###########################  LOWER MANTLE GUESSES  ###########################
	else if (model == "LOWM")
	{
		//# Activation volume 						m^3 mol^-1
		if (par == "V") val = 1e-6;
	}


return val;
}
