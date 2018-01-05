//
//  dcDataFrame.h
//  LocalSTI
//
//  Created by David CHAMPREDON on 2014-08-18. asasas
//

#ifndef __LocalSTI__dcDataFrame__
#define __LocalSTI__dcDataFrame__

#include <iostream>
#include "dcTools.h"
#include "dcMatrix.h"

#endif /* defined(__LocalSTI__dcDataFrame__) */


using namespace std;

class dcDataFrame
{

	/// Simple Data frame (similar to R)
	///
	/// For example:
	///				| _colname[0]	| _colname[1]
	///	--------------------------------------------
	/// _rowname[0]	|	1.23		|	4.56
	/// _rowname[1]	|	5.67		|	6.78
	/// _rowname[2]	|	8.79		|	3.84
	///	--------------------------------------------
	///
	
	vector<string>	_rowname;
	vector<string>	_colname;
	dcMatrix			_value;



public:
    
	
	// ======================
    // ==== CONSTRUCTORS ====
    // ======================
	
	
	dcDataFrame(){}
	
	
	dcDataFrame(vector<string> rowname, dcMatrix M, vector<string> colName)
	{
		if(rowname.size()!=M.getNbRows())
		{
			cout << "ERROR dcDataFrame constructor:";
			cout << "'_rowname' vector size and dcMatrix nb rows 'value' do not match"<<endl;
			exit(1);
		}
		
		if(colName.size()!=M.getNbCols())
		{
			cout << "ERROR dcDataFrame constructor:";
			cout << "'_colname' (headers) vector size and dcMatrix cols 'value' do not match"<<endl;
			exit(1);
		}
		
		
		_rowname = rowname;
		_value = M;
		_colname = colName;
		
	}
	
	
	dcDataFrame(vector<string> rowname, dcMatrix M)
	{
		if(rowname.size()!=M.getNbRows())
		{
			cout << "ERROR dcDataFrame constructor:";
			cout << "'_rowname' vector size and dcMatrix nb rows 'value' do not match"<<endl;
			exit(1);
		}
		
		_rowname = rowname;
		_value = M;
		
		vector<string> tmp(0);
		for (int j=0;j<_value.getNbCols(); j++)
		{
			tmp.push_back("C"+int2string(j));
		}

		_colname = tmp;
	}
	
	
	dcDataFrame(dcMatrix M)
	{
		/// Construct a dcDataFrame from a dcMatrix
		/// (row and column names are given defualt values)
		
		unsigned long ncol = M.getNbCols();
		unsigned long nrow = M.getNbRows();
		stopif(nrow==0 || ncol==0, "Cannot construct  a dcDataFrame from empty matrix");
		
		_value = M;
		
		// Default column names
		vector<string> tmp(0);
		for (int j=0;j<ncol; j++)
		{
			tmp.push_back("C"+int2string(j));
		}
		_colname = tmp;
		
		
		// Default row names
		tmp.clear();
		for (int j=0;j<nrow; j++)
		{
			tmp.push_back("R"+int2string(j));
		}
		_rowname = tmp;
	}

	
	// ==== MANIPULATIONS ====
	
	
	void addrow(string varname, vector<double> values);
	void addrow(string varname,double value);
	
	
	// ==== SET FUNCTIONS ====
	
	void set_rowname(vector<string> x){_rowname=x;}
	void set_colname(vector<string> x){_colname=x;}
	
	
	// ==== RETRIEVE DATA ====
	
	vector<string>	get_rowname() {return _rowname;}
	vector<string>	get_colname() {return _colname;}
	dcMatrix			get_value() {return _value;}
	
	
	double	getValue(string varname, string valuename);
	double	getValue(string varname, unsigned long j);
	double	getValue(string varname);
	
	
	// ==== FILE MANIPULATION ====
	void	saveToCSV(string filename,bool headers);
	
	// ==== DISPLAY ====

	void display();
	
	
};

// =================================
// ===== OUTSIDE CLASS ======
// =================================

void thetest(dcDataFrame& d);



//void to_paste_in_main_cpp()
//{
//	dcDataFrame D;
//	
//	D.addrow("test", 0.2);
//	
//	//D.display();
//	
//	dcDataFrame DD;
//	vector<double> xx,yy;
//	xx.push_back(1.234);
//	xx.push_back(5.678);
//	
//	std::vector<double>::iterator it;
//	it = xx.begin();
//	it = xx.insert ( it , 99.9 );
//	displayVector(xx);
//
//}
//
