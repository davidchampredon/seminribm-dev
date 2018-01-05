//
//  dcTools.h
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.

// Update : 2013-07-18

#ifndef dcTools_h
#define dcTools_h

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <sys/time.h>

// My Libraries
#include "dcMatrix.h"
#include "RV.h"

using namespace std;


void stopif(bool condition, string error_msg,
			int error_code=1, const char ff[]=__FUNCTION__);


int factorial(int i);
long combination(int n, int k);

void	coutline(unsigned int n);


// ===================================================
// ================ FILES MANIPULATION ===============
// ===================================================




typedef vector <double> record_t;
typedef vector <record_t> data_t;

istream& operator >> ( istream& ins, record_t& record );
istream& operator >> ( istream& ins, data_t& data );

// ===================================================
// ================ VECTOR OPERATIONS ================
// ===================================================


// Note: Templates MUST be declared in the header file (not in .cpp)

template <class T> void displayVector(vector<T> v)
{
	cout << endl<< "(size="<<v.size()<<")"<<endl<<"[";
	for (int i=0; i<v.size()-1; i++) 
	{
		cout << v[i] << "; ";
		if ((i+1)%10==0) cout<<endl;
	}
	cout << v[v.size()-1];
	cout<< "]" << endl;	
}


template <class T> vector<T> substractVector(vector<T> a, vector<T> b)
{
	assert(a.size()==b.size());
	
	vector<T> res(a.size());
	
	for (unsigned int i=0; i<a.size(); i++)
	{
		res[i] = a[i]-b[i];
	}
	
	return res;
}

double	maxElementVector(vector<double> x);
double	minElementVector(vector<double> x);
int		argminElementVector(vector<double> x);


template <class T> T sumElements(vector<T> x) 
{
	T s = 0;
	for (int i=0; i<x.size(); i++) 
	{
		s += x[i];
	}
	return s;
}

template <class T> T extractElementRandom(vector<T> x)
{
	// Extracts an element randomly
	int rndPos = uniformInt(0, (int)x.size()-1);
	
	return x[rndPos];
}


template <class T> vector<T> deleteElement(vector<T> x, int positionElementToDelete) 
{
	// Rather use "erase()" method of std vector library
	
	vector<T> y;
	for (int i=0; i<x.size(); i++) 
	{
		if (i!=positionElementToDelete)
			y.push_back(x[i]);
	}
	return y;
}


/// CHECK IF A GIVEN VALUE IS AN ELEMENT OF A VECTOR
/// (if vector empty, returns "false" anyway)
template <class T> bool isElementPresent(vector<T> x, T elemValue)
{
	bool isPresent = false;
	if (x.size()>0) 
	{
		for (unsigned long i=0; i<x.size(); i++)
		{
			if (x[i]==elemValue)
			{
				isPresent=true;
				break;
			}
		}
	}
	return isPresent;
}


template <class T> unsigned long findIndexElement(vector<T> x, T elemValue) 
{
	unsigned long s = x.size(); // initialize to a value that crashes if somethin goes wrong.
	bool isPresent = false;
	
	for (unsigned long i=0; i<x.size(); i++) 
	{
		if (x[i]==elemValue)
		{
			s=i;
			isPresent=true;
			break;
		}
	}
	
	if (!isPresent) 
	{
		cout << endl << "ERROR [findIndexElement]: element <"<< elemValue;
		cout << "> not found in this vector:";
		displayVector(x);
		exit(1);
	}
	
	return s;
}

/// Erase the first element of a vector that has element value = 'valueElementToDelete'
template <class T> vector<T> popElementValue(vector<T> x, T valueElementToDelete) 
{
    vector<T> y=x;
	unsigned long ii = findIndexElement(x, valueElementToDelete);
    y.erase(y.begin()+ii);
	return y;
}


template <class T> vector<T> popAllElementValue(vector<T> x, T valueElementToDelete)
{
	vector<T> y;
	for(int i=0;i<x.size();i++)
		if(x[i]!=valueElementToDelete)
			y.push_back(x[i]);
	return y;
}



// ===================================================
// ==================== CONVERSION ===================
// ===================================================

string int2string(int i);




#endif
