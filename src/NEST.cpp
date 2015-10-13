#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <map>
using namespace std;

#define check(flag,msg) if(flag){cerr << msg << endl;exit(1);}


// add expression contribution of a to b
void neighbor_expression_add(
		const map<string, double*> &expression_map, const size_t N,
		const string ida, const string idb, const double w,
		map<string, double*> &result_map,	// output result map
		map<string, double> &degree_map	// number of neighbors with expression
		)
{
	size_t i;
	double *result, *data;

	map<string, double*>::const_iterator eiter;
	map<string, double*>::iterator riter;

	eiter = expression_map.find(ida);
	if(eiter == expression_map.end()) return;
	data = eiter->second;

	riter = result_map.find(idb);

	if(riter == result_map.end())
	{
		result = new double[N];
		memset(result, 0, N*sizeof(double));
		riter = result_map.insert(pair<string, double*>(idb, result)).first;
		degree_map.insert(pair<string,size_t>(idb,0));
	}

	degree_map[idb] += w;

	result = riter->second;
	for(i=0;i<N;i++) result[i] += w*data[i];
}

// basic function of calculating neighbor expression from network
void neighbor_expression(
		const map<string, double*> &expression_map, const size_t N,
		const string &network,
		map<string, double*> &result_map,	// output result map
		map<string, double> &degree_map	// number of neighbors with expression
		)
{
	double w;

	string line, ida, idb;
	istringstream iss;

	ifstream fin(network.c_str());

	check(fin.fail(), "Error: Cannot open file \"" << network << "\".")

	while(getline(fin, line, '\n'))
	{
		iss.str(line);

		getline(iss, ida, '\t');
		getline(iss, idb, '\t');

		check(iss.fail(), "Error: Fail reading identifier interaction line: \"" << line << "\".")

		if(ida == idb)
		{
			iss.clear();
			continue;
		}

		if(iss.eof()) w = 1;
		else{
			iss >> w;
			check(iss.fail(), "Error: Fail reading weight in interaction line: \"" << line << "\".")
		}

		iss.clear();

		neighbor_expression_add(expression_map, N, ida, idb, w, result_map, degree_map);
		neighbor_expression_add(expression_map, N, idb, ida, w, result_map, degree_map);
	}

	fin.close();
}


void dump_map(const map<string, double*> &result_map, const size_t N, const string &header, const string &output)
{
	size_t i;
	map<string, double*>::const_iterator miter;

	ofstream fout(output.c_str());
	fout << header << '\n';

	for(miter=result_map.begin(); miter != result_map.end(); miter++)
	{
		fout << miter->first;
		for(i=0; i<N ; i++) fout << '\t' << miter->second[i];
		fout << '\n';
	}

	fout.close();
}

int main(int argc, char *argv[])
{
	// expression map
	map<string, double*> expression_map, result_map;
	map<string, double*>::iterator miter;

	map<string, double> degree_map;
	map<string, double>::iterator diter;

	string line, header, type, network, matrix, output;
	istringstream iss;
	ostringstream oss;

	size_t i, N, parseCnt = (argc-1)/2;

	double *data;

	if ( argc < 7 )
	{
		if(argc == 2 && string(argv[1]) == "-help")
		{
			cout << "\nNEST: Network Essentiality Scoring Tool\n" << endl;

			cout << "Usage: NEST -i network -m matrix -o output\n"<<endl;

			cout << "\nReport bugs to Peng Jiang (peng.jiang.software@gmail.com)\n" <<endl;
			exit(0);
		}else{
			cerr << "Insufficient number of arguments, do \"NEST -help\" for help."<<endl;
			exit(1);
		}
	}


	// read in all parameters
	for(i=0;i<parseCnt;i++)
	{
		type = argv[2*i+1];
		line = argv[2*i+2];

		if(type == "-i"){
			network = line;

		}else if(type == "-m"){
			matrix = line;

		}else if (type == "-o"){
			output = line;

		}else if (type == "-help"){
			cerr << "Please don't use \"-help\" as parameter input." << endl;
			exit(1);

		}else{
			cerr << "Cannot recognize parameter \""<< type << "\"." << endl;
			exit(1);
		}
	}

	check(network.empty(), "Cannot find input network.")
	check(matrix.empty(), "Cannot find input matrix.")
	check(output.empty(), "Cannot find output name.")

	ifstream fin;

	// load expression file
	fin.open(matrix.c_str());

	check(fin.fail(), "Cannot open file \"" << matrix << "\".")

	// read header and count number of samples
	getline(fin, header, '\n');

	iss.str(header);
	for(N=0 ; getline(iss, line, '\t') ; N++);
	iss.clear();

	while(getline(fin, line, '\n'))
	{
		iss.str(line);

		getline(iss, line, '\t');

		miter = expression_map.find(line);

		check(miter != expression_map.end(), "Error: duplicate row ID \"" << line << "\".")

		data = expression_map[line] = new double[N];

		// read in data
		for(i=0 ; i<N ; i++)
		{
			iss >> data[i];
			check(iss.fail(), "Error: Fail reading number in expression line: \"" << line << "\".")
		}

		iss.clear();
	}
	fin.close();
	fin.clear();


	// calculate neighbor expression for original network
	neighbor_expression(expression_map, N, network, result_map, degree_map);

	check(degree_map.empty(), "Empty result. Check your IDs are matched between input network and matrix.")


	dump_map(result_map, N, header, output);

	ofstream fout((output + ".degree").c_str());
	fout << "Degree\n";
	for(diter=degree_map.begin(); diter != degree_map.end(); diter++)
	{
		fout << diter->first << '\t' << diter->second << '\n';
	}
	fout.close();

	for(miter = expression_map.begin(); miter != expression_map.end(); miter++)
		delete[] miter->second;

	for(miter = result_map.begin(); miter != result_map.end(); miter++)
		delete[] miter->second;

	return 0;
}
