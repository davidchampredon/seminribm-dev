//
//  dcMatrix.h
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//

#ifndef dcMatrix_h
#define dcMatrix_h

#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include <cstdlib>


using namespace std;

class dcMatrix
{
public:
	unsigned long nbRows;	//Dimensions
	unsigned long nbCols;
    
	vector<double> val;	//Value for each element
    
	double & operator()(unsigned long,unsigned long);
	double & operator[](unsigned long);
    
    unsigned long getNbRows()     {return nbRows;}
    unsigned long getNbCols()     {return nbCols;}
    
    
    // Constructors
	
	dcMatrix(){}
	
	dcMatrix(unsigned long h,unsigned long l) 
	{
		nbRows=h; nbCols=l; 
		vector<double> tmp(l*h,0.0);
		val=tmp;
	}
	
	dcMatrix(unsigned long n) 
	{
		nbRows=n; nbCols=n;
		vector<double> tmp(n*n, 0.0);
		val=tmp;
	}
	
    dcMatrix(string pathFile);
	
	dcMatrix(vector<double> v);	// dcMatrix (n,1) from a single vector(n)
	
	
    
    void    resize(unsigned long n)  
	{
		nbRows=n;nbCols=n; 
		val.resize(n*n);
	}
    
	void    resize(unsigned long nRows, unsigned long nCols)  
	{
		nbRows=nRows; nbCols=nCols; 
		this->val.resize(nRows*nCols);
	}
	
	void	clear() {val.clear();nbCols=0; nbRows=0;}
    
	void    display();

    
    // Files operations
	void    FromFile(const char*);
    void    FromFile(string);
	void	FromFile_Rows(string fileName, unsigned long nrow);
	
    void    WriteToFile(string);
	void    WriteToFileCSV(string);
	void    WriteToFileCSV(string fileName, vector<string> header);
	
	
    // Operations on embedded vectors
    vector<double>  extractColumn(unsigned long j_col);
	vector<double>	extractRow(unsigned long i_row);
	
	void            addRowVector(vector<double> v);
	void            addRowVector(vector<unsigned long> v);
	
    void            addColVector(vector<double> v);
    
	void			removeRow(unsigned long i_row);	// removes row 'i_row' and resize dcMatrix
	void			removeCol(unsigned long j_col);	// removes column 'j_col' and resize dcMatrix

	
	
	// Extract the row #i of the matrix
	// "i" is calculated such that it is
	// the smallest element of column "j_col"
  	vector<double>	extractRow_cond_minElement(unsigned long j_col);
	
    // Operations on elements
    
    void    setAllValues(double value);
	void	setValueFromMatrix(dcMatrix M);
	
	double  sumAllElements();
	double  sumLine(unsigned long i);		// sum all elements of line #i
	double  sumColumn(unsigned long j);	// sum all elements of line #i
	
	// conditional sum 
	double	sumColumnIf(unsigned long colToSum, unsigned long colToTest,
						double lowerBound, double upperBound);

	// counts nb of elements which are lower<element<upper  
	unsigned long		countColumnIf(unsigned long colToTest,
						  double lowerBound, double upperBound);
    
    void    setColumnValues(unsigned long colNb, vector<double> v);
	void	setRowValues(unsigned long rowNb_start0, vector<double> v);
	void	setRowValues(unsigned long rowNb_start0, vector<unsigned long> v);
	
	unsigned long countNonZeroElements();
	unsigned long countNonZeroElements_line(unsigned long i);
	
    dcMatrix  transpose();
    
    bool    isSymetric();
    
    double  getMinimumValue();
    double  getMaximumValue();
	
};

dcMatrix operator + (dcMatrix &A,dcMatrix &B);
dcMatrix operator - (dcMatrix &A,dcMatrix &B);
dcMatrix operator * (dcMatrix &A,dcMatrix &B);
dcMatrix operator * (double a,dcMatrix &A);

dcMatrix Id(unsigned long n);

dcMatrix power(dcMatrix A,int n);

double distance_Matrix(dcMatrix A, dcMatrix B, double power);	// Euclidian distance b/w two matrices

dcMatrix rowBind(dcMatrix A, dcMatrix B);



#endif
