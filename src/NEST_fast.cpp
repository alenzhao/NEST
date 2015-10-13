#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <climits>
#include <map>
#include <cmath>
using namespace std;

typedef unsigned short int USHRT;

#define EPS 1e-10
#define check(flag,msg) if(flag){cerr << msg << endl;exit(1);}


// add expression contribution of a to b
inline void neighbor_expression_add(
		const double* const expression_map[], const size_t N,
		const USHRT inxa, const USHRT inxb,
		double* result_map[],	// output result map
		size_t degree_map[]	// number of neighbors with expression
		)
{
	const double* const data = expression_map[inxa];
	if(!data) return;

	if(degree_map) degree_map[inxb] ++;

	double *result = result_map[inxb];

	for(size_t i=0;i<N;i++) result[i] += data[i];
}


// basic function of calculating neighbor expression from network
void neighbor_expression(
		const double* const expression_map[],
		const size_t N,
		const string &network,
		double* result_map[],
		size_t degree_map[])
{
	USHRT inxa, inxb;

	ifstream fin(network.c_str());
	check(fin.fail(), "Error: Cannot open file \"" << network << "\".")

	while(fin.read((char *)(&inxa), sizeof(USHRT)) && fin.read((char *)(&inxb), sizeof(USHRT)))
	{
		if(inxa == inxb) continue;

		neighbor_expression_add(expression_map, N, inxa, inxb, result_map, degree_map);
		neighbor_expression_add(expression_map, N, inxb, inxa, result_map, degree_map);
	}

	fin.close();
}


void dump_map(
		const double* const result_map[],
		const size_t degree_map[],
		const size_t Ngene, const size_t Nsample,
		const vector<string> name_vec, const string &header,
		const string &output)
{
	size_t i, j;

	ofstream fout(output.c_str());
	fout << header << '\n';

	for(i=0;i<Ngene;i++)
	{
		if(!degree_map[i]) continue;

		fout << name_vec[i];
		for(j=0; j<Nsample ; j++) fout << '\t' << result_map[i][j];
		fout << '\n';
	}

	fout.close();
}


int main(int argc, char *argv[])
{
	vector<string> name_vec;
	map<string, USHRT> name_map;
	map<string, USHRT>::iterator miter;

	string line, header, type, network, matrix, output;
	istringstream iss;
	ostringstream oss;

	size_t i, j, Nsample, Ngene, Nbuffer, jumped, *degree_map, *p_l_arr, *p_g_arr, Nrand = 1000, parseCnt = (argc-1)/2;

	bool write_rand_result = false, flag_empty = true;

	char *pEnd;
	double v, *data, **expression_map,
		**result_map, *result_arr,
		**rand_map, *rand_arr,
		*sum_arr, *sumsq_arr;

	if ( argc < 7 )
	{
		if(argc == 2 && string(argv[1]) == "-help")
		{
			cout << "\nNEST_fast: Network Essentiality Scoring Tool with statistical inference\n" << endl;

			cout << "Usage: NEST_fast -i network -m matrix -o output [OPTIONS]\n"<<endl;
			cout << "\tnetwork: prefix of compacted binary random networks\n" << endl;

			cout << "Options:" << endl;
			cout << "\t-n\tNumber of randomizations. Default: " << Nrand << endl;
			cout << "\t-w\tWrite intermediate results: 1 (yes) or 0 (no). Default: " << (write_rand_result?"1 (yes)":"0 (no)") << endl;
			cout << "\nReport bugs to Peng Jiang (peng.jiang.software@gmail.com)\n" <<endl;
			exit(0);
		}else{
			cerr << "Insufficient number of arguments, do \"NEST_fast -help\" for help."<<endl;
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

		}else if (type == "-n"){
			Nrand = strtol(line.c_str(), &pEnd, 10);
			check(*pEnd!='\0' || pEnd == line.c_str(), "Error: \"" << line << "\" is not an integer.")

		}else if (type == "-w"){
			check(line != "0" && line != "1", "Please input 0 or 1 for flag. \"" << line << "\" is input.")
			write_rand_result = (line != "0");

		}else if (type == "-help"){
			cerr << "Please don't use \"-help\" as parameter input." << endl;
			exit(1);

		}else{
			cerr << "Cannot recognize parameter \""<< type << "\"." << endl;
			exit(1);
		}
	}

	check(network.empty(), "Error: cannot find input network.")
	check(matrix.empty(), "Error: cannot find input matrix.")
	check(output.empty(), "Error: cannot find output name.")

	ifstream fin;

	// load name map

	fin.open((network + ".map").c_str());
	check(fin.fail(), "Error: Cannot open name map \"" << network << ".map\".")

	for(i=0;getline(fin, line, '\n');i++)
	{
		name_map[line] = i;
		name_vec.push_back(line);
	}

	fin.close();
	fin.clear();

	Ngene = name_vec.size();

	check(Ngene > USHRT_MAX, "Error: The index " << Ngene << " exceeds limit " << USHRT_MAX)

	// load expression file
	fin.open(matrix.c_str());

	check(fin.fail(), "Cannot open file \"" << matrix << "\".")

	// allocate pointer space for expression
	expression_map = new double*[Ngene];
	memset(expression_map, 0, Ngene*sizeof(double*));

	// read header and count number of samples
	getline(fin, header, '\n');

	iss.str(header);
	for(Nsample=0 ; getline(iss, line, '\t') ; Nsample++);
	iss.clear();

	Nbuffer = Ngene*Nsample;

	// Number of sample is know, allocate space for result
	degree_map = new size_t[Ngene];
	result_map = new double*[Ngene];
	result_arr = new double[Nbuffer];

	memset(degree_map, 0, Ngene*sizeof(size_t));
	memset(result_arr, 0, Nbuffer*sizeof(double));

	for(i=j=0; i<Ngene; i++,j+=Nsample)
	{
		result_map[i] = result_arr + j;
	}

	jumped = 0;

	while(getline(fin, line, '\n'))
	{
		iss.str(line);

		getline(iss, line, '\t');

		miter = name_map.find(line);

		if(miter == name_map.end())
		{
			jumped ++;
			iss.clear();
			continue;
		}

		check(expression_map[miter->second], "Error: duplicate row ID \"" << line << "\".")

		data = expression_map[miter->second] = new double[Nsample];

		// read in data
		for(i=0 ; i<Nsample ; i++)
		{
			iss >> data[i];
			check(iss.fail(), "Error: Fail reading number in expression line: \"" << line << "\".")
		}

		iss.clear();
	}
	fin.close();
	fin.clear();

	if(jumped) cout << jumped << " genes in matrix jumped (not included in network)." << endl;


	// calculate neighbor expression for original network
	neighbor_expression(expression_map, Nsample, network, result_map, degree_map);

	dump_map(result_map, degree_map, Ngene, Nsample, name_vec, header, output);

	ofstream fout((output + ".degree").c_str());
	fout << "Degree\n";
	for(i=0; i<Ngene; i++)
		if(degree_map[i])
		{
			fout << name_vec[i] << '\t' << degree_map[i] << '\n';
			flag_empty = false;
		}
	fout.close();

	check(flag_empty, "Empty result. Check your IDs are matched between input network and matrix.")

	if(Nrand > 0)
	{
		// create space for randomization result

		rand_map = new double*[Ngene];
		rand_arr = new double[Nbuffer];
		sum_arr = new double[Nbuffer];
		sumsq_arr = new double[Nbuffer];

		p_l_arr = new size_t[Nbuffer];
		p_g_arr = new size_t[Nbuffer];

		memset(p_l_arr, 0, Nbuffer * sizeof(size_t));
		memset(p_g_arr, 0, Nbuffer * sizeof(size_t));
		memset(sum_arr, 0, Nbuffer * sizeof(double));
		memset(sumsq_arr, 0, Nbuffer * sizeof(double));

		// hook up result ID to random maps
		for(i=j=0; i< Ngene; i++, j+=Nsample) rand_map[i] = rand_arr + j;

		cout << "process random networks" << endl;

		for(i=0;i<Nrand;i++)
		{
			cout << i << endl;
			oss << network << '.' << i;

			memset(rand_arr, 0, Nbuffer * sizeof(double));
			neighbor_expression(expression_map, Nsample, oss.str(), rand_map, NULL);
			oss.str(""); oss.clear();

			if(write_rand_result)
			{
				oss << output << '.' << i;
				dump_map(rand_map, degree_map, Ngene, Nsample, name_vec, header, oss.str());
				oss.str(""); oss.clear();
			}

			for(j=0; j<Nbuffer; j++)
			{
				v = rand_arr[j];

				if(v >= result_arr[j]) p_g_arr[j] ++;
				if(v <= result_arr[j]) p_l_arr[j] ++;

				sum_arr[j] += v;
				sumsq_arr[j] += v*v;
			}
		}

		// calculate Z-score and P-value
		for(i=0;i<Nbuffer;i++)
		{
			v = sum_arr[i] = sum_arr[i]/Nrand;
			v = sumsq_arr[i]/Nrand - v * v;

			// use sum_arr as Z, sumsq_arr as P-value
			if(fabs(v) < EPS){
				sum_arr[i] = 0;
				sumsq_arr[i] = 1;
			}else{
				check(v < 0, "Error: variation is smaller than zero. Please report bug to peng.jiang.software@gmail.com")

				  v = sum_arr[i] = (result_arr[i] - sum_arr[i])/sqrt(v);

				if(v>0)
					sumsq_arr[i] = (double)p_g_arr[i]/Nrand;
				else
					sumsq_arr[i] = (double)p_l_arr[i]/Nrand;
			}
		}

		for(i=j=0; i< Ngene; i++, j+=Nsample) rand_map[i] = sum_arr + j;
		dump_map(rand_map, degree_map, Ngene, Nsample, name_vec, header, output + ".Zscore");

		for(i=j=0; i< Ngene; i++, j+=Nsample) rand_map[i] = sumsq_arr + j;
		dump_map(rand_map, degree_map, Ngene, Nsample, name_vec, header, output + ".Pvalue");

		delete[] rand_map; delete[] rand_arr;
		delete[] p_l_arr; delete[] p_g_arr;
	}

	for(i=0;i<Ngene;i++)
		if(expression_map[i])
			delete[] expression_map[i];

	delete[] expression_map;
	delete[] degree_map;
	delete[] result_map; delete[] result_arr;

	return 0;
}
