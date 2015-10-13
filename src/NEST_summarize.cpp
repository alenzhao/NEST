#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
using namespace std;

#define check(flag,msg) if(flag){cerr << msg << endl;exit(1);}


struct rank_node
{
	double v; size_t i;
	bool operator < (const rank_node &r) const {return v < r.v;}
};

void Benjamini_Hochberg(
		const double pvalue[], const size_t N,	// input
		double FDR[],	// output
		rank_node sort_arr[]		// workspace
	)
{
	size_t i;
	double qvalue = 1;

	for(i=0;i<N;i++)
	{
		sort_arr[i].v = pvalue[i];
		sort_arr[i].i = i;
	}

	sort(sort_arr, sort_arr+N);

	for(i=0;i<N;i++)
	{
		qvalue = min<double>(qvalue, sort_arr[N-1-i].v*N/(N-i));
		FDR[sort_arr[N-1-i].i] = qvalue;
	}
}


int main(int argc, char *argv[])
{
	size_t i, j, Ncol, Nrow, Nbuf;
	string line;
	istringstream iss;
	vector<string> colnames, rownames;
	vector<double> *vvec = new vector<double>();

	double v, *pvalue, *qvalue;
	rank_node *sort_arr;

	ifstream fin, fin_degree, fin_zscore;
	ofstream fout;

	if(argc != 2)
	{
		cerr << "NEST_summarize output_prefix" << endl;
		return 1;
	}

	string output_prefix = argv[1];

	fin.open((output_prefix + ".Pvalue").c_str());
	check(fin.fail(), "Cannot open file \"" << output_prefix << ".Pvalue\"")

	// read header and count number of columns
	getline(fin, line, '\n');

	iss.str(line);
	while(getline(iss, line, '\t')) colnames.push_back(line);
	iss.clear();

	Ncol = colnames.size();


	while(getline(fin, line, '\n'))
	{
		iss.str(line);
		getline(iss, line, '\t');
		rownames.push_back(line);

		for(i=0;i<Ncol;i++)
		{
			iss >> v;
			vvec->push_back(v);
		}

		iss.clear();
	}

	fin.close();
	fin.clear();

	Nrow = rownames.size();
	Nbuf = Ncol * Nrow;

	pvalue = new double[Nbuf];

	// transpose all P-values
	for(i=0;i<Nrow;i++)
		for(j=0;j<Ncol;j++)
			pvalue[j*Nrow + i] = vvec->at(i*Ncol+j);

	delete vvec;
	qvalue = new double[Nbuf];

	// workspace for Benjamini Hochberg correction
	sort_arr = new rank_node[Nrow];

	for(i=0;i<Ncol;i++) Benjamini_Hochberg(pvalue + i*Nrow, Nrow, qvalue + i*Nrow, sort_arr);


	fout.open((output_prefix + ".summary").c_str());

	fout << "ID\tDegree";

	for(i=0;i<Ncol;i++)
		fout << '\t' << colnames[i]
		     << '\t' << "Zscore"
		     << '\t' << "Pvalue"
		     << '\t' << "FDR";

	fout << '\n';

	fin.open(output_prefix.c_str());
	check(fin.fail(), "Cannot open \"" << output_prefix << "\"")
	getline(fin,line,'\n');

	fin_degree.open((output_prefix + ".degree").c_str());
	check(fin_degree.fail(), "Cannot open \"" << output_prefix << ".degree\"")
	getline(fin_degree,line,'\n');

	fin_zscore.open((output_prefix + ".Zscore").c_str());
	check(fin_zscore.fail(), "Cannot open \"" << output_prefix << ".Zscore\"")
	getline(fin_zscore,line,'\n');

	for(i=0;i<Nrow;i++)
	{
		// check the row names are consistent
		getline(fin, line, '\t');
		check(line != rownames[i], "Mismatch row names between \"" << line << "\" and \"" << rownames[i] << "\" in score file.");

		getline(fin_degree, line, '\t');
		check(line != rownames[i], "Mismatch row names between \"" << line << "\" and \"" << rownames[i] << "\" in degree file.");

		getline(fin_zscore, line, '\t');
		check(line != rownames[i], "Mismatch row names between \"" << line << "\" and \"" << rownames[i] << "\" in Z-score file.");

		fin_degree >> v;
		fout << rownames[i] << '\t' << v;

		for(j=0;j<Ncol;j++)
		{
			fin >> v; fout << '\t' << v;

			fin_zscore >> v; fout << '\t' << v;

			fout << '\t' << pvalue[j*Nrow+i] << '\t' << qvalue[j*Nrow+i];
		}

		// jump the rear symbol
		getline(fin, line, '\n');
		getline(fin_degree, line, '\n');
		getline(fin_zscore, line, '\n');

		fout << '\n';
	}

	fin_degree.close();
	fin_zscore.close();
	fin.close();
	fout.close();

	delete[] sort_arr;
	delete[] pvalue;
	delete[] qvalue;

	return 0;
}
