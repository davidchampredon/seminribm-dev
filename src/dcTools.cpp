//
//  dcTools.cpp
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//

#include "dcTools.h"



void stopif(bool condition, string error_msg,
			int error_code, const char ff[])
{
	if (condition)
	{
		cerr << endl << " *=*=*=*=*=*=* ERROR *=*=*=*=*=*=* " << endl<<endl;
		cerr<<"In function: "<<string(ff)<<endl<<endl;
		cerr << error_msg <<endl;
		cerr <<	endl <<	" *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* " << endl;
		exit(error_code);
	}
}

int factorial(int i) 
{
    if (i <= 1)
        return i;
    return (i * factorial(i - 1));
}


long combination(int n, int k)
{
	assert(n>=k);
	return factorial(n)/factorial(k)/factorial(n-k);
}



void coutline(unsigned int n)
 {cout<<endl; for (int i=0;i<n;i++) cout<<"_"; cout<<endl;}


// ===================================================
// ================ FILES MANIPULATION ===============
// ===================================================



typedef vector <double> record_t;
typedef vector <record_t> data_t;

typedef vector <string> record_t_string;
typedef vector <record_t_string> data_t_string;


//-----------------------------------------------------------------------------
// Let's overload the stream input operator to read a list of CSV fields (which a CSV record).
// Remember, a record is a list of doubles separated by commas ','.

istream& operator >> ( istream& ins, record_t& record )
{
	
	record.clear();	// make sure that the returned record contains only the stuff we read now
	
	string line;
	getline( ins, line );	// read the entire line into a string (a CSV record is terminated by a newline)
    
	stringstream ss( line );		// now we'll use a stringstream to separate the fields out of the line
	string field;
    
    //double tmp; int n;
    
	while (getline( ss, field, ',' ))
    {
		// for each field we wish to convert it to a double
		// (since we require that the CSV contains nothing but floating-point values)
		stringstream fs( field );
		double f = 0.0;  // (default value is 0.0)
		fs >> f;
		
		record.push_back( f );		// add the newly-converted field to the end of the record
        
		//n=record.size();       // FIXME: to delete
		//tmp = record[n-1];
    }
	
	// Now we have read a single line, converted into a list of fields, converted the fields
	// from strings to doubles, and stored the results in the argument record, so
	// we just return the argument stream as required for this kind of input overload function.
	
	return ins;
}

istream& operator >> ( istream& ins, record_t_string& record )
{
	
	record.clear();	// make sure that the returned record contains only the stuff we read now
	
	string line;
	getline( ins, line );	// read the entire line into a string (a CSV record is terminated by a newline)
    
	stringstream ss( line );		// now we'll use a stringstream to separate the fields out of the line
	string field;
    
	while (getline( ss, field, ',' ))
    {
		// for each field we wish to convert it to a double
		// (since we require that the CSV contains nothing but floating-point values)
		stringstream fs( field );
		string f = "";  // (default value is 0.0)
		fs >> f;
		
		record.push_back( f );		// add the newly-converted field to the end of the record
        
		//n=record.size();       // FIXME: to delete
		//tmp = record[n-1];
    }
	
	// Now we have read a single line, converted into a list of fields, converted the fields
	// from strings to doubles, and stored the results in the argument record, so
	// we just return the argument stream as required for this kind of input overload function.
	
	return ins;
}

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.
// This time it is a little easier, just because we only need to worry about reading
// records, and not fields.

istream& operator >> ( istream& ins, data_t& data )
{
	data.clear();
	
	// For every record we can read from the file, append it to our resulting data
	record_t record;
	while (ins >> record)
    {
		data.push_back( record );
    }
	
	return ins;  
}

istream& operator >> ( istream& ins, data_t_string& data )
{
	data.clear();
	
	// For every record we can read from the file, append it to our resulting data
	record_t_string record;
	while (ins >> record)
    {
		data.push_back( record );
    }
	
	return ins;
}


// ===================================================
// ================ VECTOR OPERATIONS ================
// ===================================================




double maxElementVector(vector<double> x)
{
    double max=-999999999;
    
    for (int i=0; i<x.size(); i++) {
        if (max < x[i]) {
            max = x[i];
        }
    }
    return max;
}

double minElementVector(vector<double> x)
{
    double min = 999999999;
    
    for (int i=0; i<x.size(); i++) 
	{
        if (min > x[i]) 
		{
            min = x[i];
        }
    }
    return min;
}

int argminElementVector(vector<double> x)
{
    double min = 999999999;
	int argmin = 0;
    
    for (int i=0; i<x.size(); i++) 
	{
        if (min > x[i]) 
		{
            min = x[i];
			argmin = i;
        }
    }
    return argmin;
}




// ===================================================
// ==================== CONVERSION ===================
// ===================================================



string int2string(int i)
{
	string res;					// string which will contain the result
	ostringstream convert;		// stream used for the conversion
	convert << i;				// insert the textual representation of 'Number' in the characters in the stream
	res = convert.str();		// set 'Result' to the contents of the stream
	
	return res;	
}








