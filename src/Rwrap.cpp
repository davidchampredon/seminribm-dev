//
//  Rwrap.cpp
//  naf
//
//  Created by David CHAMPREDON on 2018-01-03.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <Rcpp.h>

using namespace Rcpp;

#include "simulator.h"
#include "individual.h"
#include "dcDataFrame.h"
#include "globalVar.h"



// [[Rcpp::export]]
List cpp_seminribm_run(double horizon,
				   unsigned long popSize,
				   double R0  ,
				   double latent_mean ,
				   double infectious_mean,
				   int nE,
				   int nI,
				   unsigned long initInfectious,
				   bool calc_WIW_Re,
				   bool doExact,
				   double timeStepTauLeap,
				   unsigned int rnd_seed){
    
    try{
        cout << " SEmInR simulation launched... " << endl;
		
        // =========================
        // === Call C++ function ===
        // =========================
        
        _RANDOM_GENERATOR.seed(rnd_seed);
		
        // Derive other variables
        double sigma0    = 1.0/latent_mean;
        double gamma0    = 1.0/infectious_mean;
        double beta      = R0 * gamma0;
		
        vector<double> sigma(nE);
        vector<double> gamma(nI);
        for (int i=0; i<nE; i++) sigma[i]=sigma0*nE;
        for (int i=0; i<nI; i++) gamma[i]=gamma0*nI;
        
        // Simulation
        simulator SIM(beta, sigma, gamma, popSize, nE, nI);
        SIM.run(horizon, initInfectious, calc_WIW_Re, doExact, timeStepTauLeap);
		
		// Retrieve outputs:
        vector<double> times		= SIM.get_time();
        vector<unsigned long> ns 	= SIM.get_nS();
		vector<unsigned long> nr 	= SIM.get_nR();
		vector<unsigned long> prev 	= SIM.get_prevalence();
		vector<double> acq_t 		= SIM.get_timeDiseaseAcquisition();
		vector<double> acqTransm_t 	= SIM.get_timeDiseaseAcquisition_transm();
		vector<double> b 			= SIM.get_GIbck();
		vector<vector<double> > f 	= SIM.get_GIfwd();
		vector<vector<double> > reff = SIM.get_Reff();
		
        Rcpp::List empty_list;
        
        return List::create(Named("times") = times,
							Named("S")     = ns,
							Named("prev")  = prev,
							Named("R")     = nr,
							Named("acq_times")   = acq_t,
							Named("acqTransm_t") = acqTransm_t,
							Named("GI_bck") = b,
							Named("GI_fwd") = f,
							Named("Reff") = reff
                            );
    }
    
    catch (...){
        ::Rf_error(">>>> C++ exception (unknown reason) <<<<");
        return NULL;
    }
}


