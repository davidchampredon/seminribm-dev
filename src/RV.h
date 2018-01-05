/*
 *  RV.h
 *  Gillespie
 *
 *  Created by David Champredon  on 13-04-19.
 *
 */

#ifndef RV_H
#define RV_H

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <random>

using namespace std;


// === Random number generator seed ====

void force_seed_reset();
void force_seed_reset(unsigned int);


// === Continuous support ====

double	uniform(int seed);				// Uniform in [0;1]
double	uniform01();

double	expo(double lambda);
double	normal(double mean, double stddev);

double	gamma(double shape, double scale);
double	beta(double a, double b);


// === Discrete support ====

int		uniformInt(int min, int max); // Uniform INTEGER in [min;max]
int		geometric(double p);
int		binom(double p, int N);	// Binomial(p,N) N trials with probability of success p
unsigned long	poisson(double expectedValue);

vector<unsigned int> multinomial(unsigned int N, 
								 vector<double> proba); // Multinomial (N,p) ; returns vector of size proba.size()



#endif
