#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <map>
#include <list>
using namespace std;

#define check(flag,msg) if(flag){cerr << msg << endl;exit(1);}


bool abs_comp (const double a, const double b) { return fabs(a) > fabs(b); }

double merge(vector<double> &vec, const string &method)
{
	if(method == "mean"){
		return accumulate(vec.begin(), vec.end(), 0.0)/vec.size();

	}else if(method == "median")
	{
		size_t mid = vec.size()/2;

		sort(vec.begin(), vec.end());

		if(vec.size()%2==0)
			return (vec[mid] + vec[mid-1])/2;
		else
			return vec[mid];

	}else if(method == "second")
	{
		sort(vec.begin(), vec.end(), abs_comp);
		return vec[1];

	}else if(method == "top3")
	{
		sort(vec.begin(), vec.end(), abs_comp);
		return (vec[0] + vec[1] + vec[2])/3;

	}else{
		cerr << "Cannot recognize merge method \"" << method << "\"" << endl;
		exit(1);
	}

	return 0;
}

int main(int argc, char *argv[])
{
	size_t i, N, parseCnt = (argc-1)/2, min_count = 3;

	string line, type, input, output, merge_method = "median";
	double v;

	vector<string> header;
	vector<string>::iterator hiter;
	vector<double> vec, vec_merged;

	map<string, list<double>* > merge_map;
	map<string, list<double>* >::iterator miter;

	list<double>* lptr;
	list<double>::iterator liter;

	ifstream fin;
	istringstream iss;


	if ( argc < 5 )
	{
		if(argc == 2 && string(argv[1]) == "-help")
		{
			cout << "\nmerge_matrix: merge matrix with multiple values for the same gene\n" << endl;
			cout << "Usage: merge_matrix -i input -o output [Options]\n"<<endl;

			cout << "Options:" << endl;
			cout << "\t-m\tMerge method. Default: " << merge_method << endl;
			cout << "\t\tmedian: median of all values" << endl;
			cout << "\t\tmean: mean of all values" << endl;
			cout << "\t\ttop3: mean of top three values by absolute" << endl;
			cout << "\t\tsecond: second best values by absolute" << endl;
			cout << endl;
			cout << "\t-c\tMinimum number threshold. Default: " << min_count << endl;

			cout << "\nReport bugs to Peng Jiang (peng.jiang.software@gmail.com)\n" <<endl;
			exit(0);
		}else{
			cerr << "Insufficient number of arguments, do \"merge_matrix -help\" for help."<<endl;
			exit(1);
		}
	}

	// read in all parameters
	for(i=0;i<parseCnt;i++)
	{
		type = argv[2*i+1];
		line = argv[2*i+2];

		if(type == "-i"){
			input = line;

		}else if (type == "-o"){
			output = line;

		}else if (type == "-m"){
			merge_method = line;

		}else if (type == "-c"){
			min_count = atoi(line.c_str());

		}else if (type == "-help"){
			cerr << "Please don't use \"-help\" as parameter input." << endl;
			exit(1);
		}else{
			cerr << "Cannot recognize parameter \""<< type << "\"." << endl;
			exit(1);
		}
	}

	check(input.empty(), "Cannot find input")
	check(output.empty(), "Cannot find output")

	fin.open(input.c_str());
	check(fin.fail(), "Cannot open file \"" << input << "\"")

	// read header and count number of samples
	getline(fin, line, '\n');

	iss.str(line);

	for(N=0, getline(iss, line, '\t') ; getline(iss, line, '\t') ; N++)
	{
		header.push_back(line);
	}

	iss.clear();

	while(getline(fin, line, '\n'))
	{
		iss.str(line);
		getline(iss, line, '\t');

		miter = merge_map.find(line);

		if(miter == merge_map.end())
		{
			lptr = new list<double>[N];
			merge_map.insert(pair<string, list<double>* >(line, lptr));
		}else{
			lptr = miter->second;
		}

		// jump index field
		getline(iss, line, '\t');

		// read in data
		for(i=0 ; i<N ; i++)
		{
			iss >> v;
			check(iss.fail(), "Error: Fail reading number in expression line: \"" << line << "\".")

			lptr[i].push_back(v);
		}

		iss.clear();
	}

	fin.close();
	fin.clear();

	ofstream fout(output.c_str());

	for(i=0;i<N;i++) fout << header[i] << (i==N-1?'\n':'\t');

	for(miter = merge_map.begin(); miter != merge_map.end(); miter++)
	{
		lptr = miter->second;

		for(i=0;i<N;i++)
		{
			for(liter = lptr[i].begin(); liter != lptr[i].end(); liter++)
			{
				vec.push_back(*liter);
			}

			if(vec.size() >= min_count)
			{
				v = merge(vec, merge_method);
				vec_merged.push_back(v);
			}

			vec.clear();
		}

		if(vec_merged.size() == N)
		{
			fout << miter->first;
			for(i=0;i<N;i++) fout << '\t' << vec_merged[i];
			fout << '\n';
		}

		vec_merged.clear();

		delete[] lptr;
	}

	fout.close();

	return 0;
}
