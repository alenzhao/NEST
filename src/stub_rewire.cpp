#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
using namespace std;

#define check(flag,msg) if(flag){cerr << msg << endl;exit(1);}

typedef unsigned short int USHRT;


size_t randomize(
		const vector<size_t> &degree_vec, size_t candidate_arr[],	// input
		size_t *adjacency_list[],	// output
		size_t degree_cur[], bool flag[]	// workspace
){
	size_t i,j, inx, offset, N, Ncandidate, count, connected = 0;

	N = Ncandidate = degree_vec.size();

	memset(degree_cur, 0, N*sizeof(size_t));
	memset(flag, 0, N*sizeof(bool));

	while(Ncandidate > 1)
	{
		// take out first vertex and connect all its edges
		i = *candidate_arr;
		candidate_arr ++;
		Ncandidate--;

		// full, jump to next
		if(degree_cur[i] == degree_vec[i]) continue;

		// mark visited
		for(j=0;j<degree_cur[i];j++) flag[adjacency_list[i][j]] = true;

		// for unconnected edges of current vertex
		for(j=degree_cur[i]; j<degree_vec[i]; j++)
		{
			// circular searching for unconnected vertex
			for(offset = rand() % Ncandidate, count=0; count < Ncandidate; offset = (offset+1) % Ncandidate, count ++)
			{
				inx = candidate_arr[offset];

				if(not flag[inx] && degree_cur[inx] < degree_vec[inx]) break;
			}

			// cannot continue on this vertex
			if (count == Ncandidate) break;

			adjacency_list[i][j] = inx;
			adjacency_list[inx][degree_cur[inx]] = i;
			degree_cur[i]++;
			degree_cur[inx]++;
			connected++;

			flag[inx] = true;
		}

		// erase marks
		for(j=0;j<degree_cur[i];j++) flag[adjacency_list[i][j]] = false;
	}

	return connected;
}

void output_network(	const vector<string> &name_vec, const size_t degree_arr[],
		size_t *adjacency_list[], const string &output, const bool compact_mode)
{
	size_t i, j, inx, N = name_vec.size();
	string ida, idb;

	ofstream fout(output.c_str());

	for(i=0;i<N;i++)
	{
		ida = name_vec[i];

		for(j=0;j<degree_arr[i];j++)
		{
			inx = adjacency_list[i][j];
			idb = name_vec[inx];

			if (ida < idb)
			{
				if(compact_mode){
					fout.write((const char *)(&i), sizeof(USHRT));
					fout.write((const char *)(&inx), sizeof(USHRT));
				}else{
					fout << ida << '\t' << idb << '\n';
				}
			}
		}
	}

	fout.close();
}

int main(int argc, char *argv[])
{
	// check duplicated edges
	set<size_t> included;

	// name to index
	map<string, size_t> idmap;
	map<string, size_t>::iterator miter;

	// name vectors
	vector<string> name_vec;
	vector<size_t> degree_vec;

	string line, type, ida, idb, input, output;

	istringstream iss;
	ostringstream oss;

	double ratio, max_ratio, rthres = 0.98;

	size_t i, j, inx, Nvertex, Nedge, Nconnected, Nrand = 1000, Nretry = 10,
		**adjacency_list, *adjacency_list_edge_index, *degree_cur, *candidate_arr,
		parseCnt = (argc-1)/2;

	bool *flag, compact_mode = true;
	char *pEnd;

	if ( argc < 5 )
	{
		if(argc == 2 && string(argv[1]) == "-help")
		{
			cout << "\nstub_rewire: randomize unweighted undirected network\n" << endl;

			cout << "Usage: stub_rewire -i input -o output [OPTIONS]\n"<<endl;

			cout << "\tinput: Input network "<< endl;
			cout << "\toutput: Randomized network output prefix "<< endl;

			cout << "\nOptions:" << endl;
			cout << "\t-n\tNumber of randomizations. Default: " << Nrand << endl;
			cout << "\t-r\tRetry count. Default " << Nretry << endl;
			cout << "\t-c\tConnected threshold. Default: " << rthres << endl;
			cout << "\t-m\tCompact mode: 1 (yes) or 0 (no). Default: " << (compact_mode?"1 (yes)":"0 (no)") << endl;

			cout << "\nReport bugs to Peng Jiang (peng.jiang.software@gmail.com)\n" <<endl;
			exit(0);
		}else{
			cerr << "Insufficient number of arguments, do \"stub_rewire -help\" for help."<<endl;
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

		}else if (type == "-n"){
			// number of randomizations
			Nrand = strtol(line.c_str(), &pEnd, 10);
			check(*pEnd!='\0' || pEnd == line.c_str(), "Error: \"" << line << "\" is not an integer.")

		}else if (type == "-r"){
			// number of retries if stub-rewiring process cannot continue
			Nretry = strtol(line.c_str(), &pEnd, 10);
			check(*pEnd!='\0' || pEnd == line.c_str(), "Error: \"" << line << "\" is not an integer.")
			check(Nretry==0, "Error: retry count cannot be zero.")

		}else if (type == "-c"){
			// minimum ratio of successful connected edges in stub-rewiring
			rthres = strtod(line.c_str(), &pEnd);
			check(*pEnd!='\0' || pEnd == line.c_str(), "Error: \"" << line << "\" is not a float number.")
			check(rthres <= 0 || rthres > 1, "Error: connected threshold \"" << rthres << "\" is out of range (0,1].")

		}else if (type == "-m"){
			// whether output randomized network in compacted binary files
			check(line != "0" && line != "1", "Error: Please input 1 or 0 for compact mode. You input \"" << line << "\".")
			compact_mode = (line != "0");

		}else if (type == "-help"){
			cerr << "Error: Please don't use \"-help\" as parameter input." << endl;
			exit(1);

		}else{
			cerr << "Error: Cannot recognize parameter \""<< type << "\"." << endl;
			exit(1);
		}
	}

	check(input.empty(), "Error: cannot find input network.")
	check(output.empty(), "Error: cannot find output prefix.")

	ifstream fin(input.c_str());
	check(fin.fail(), "Cannot open file " << input)

	ofstream fout(output.c_str());

	while(getline(fin, line, '\n'))
	{
		iss.str(line);

		getline(iss, ida, '\t');
		getline(iss, idb, '\t');
		check(iss.fail(), "Fail reading identifier interaction line: \"" << line << "\".")
		iss.clear();

		if(ida == idb)
		{
			cerr << "Warning: Jump self edge on " << ida << endl;
			continue;
		}

		if(ida > idb) swap(ida, idb);

		// find ida index
		miter = idmap.find(ida);
		if(miter == idmap.end())
		{
			miter = idmap.insert(pair<string, size_t>(ida, name_vec.size())).first;
			name_vec.push_back(ida);
			degree_vec.push_back(0);
		}
		i = miter->second;

		// find idb index
		miter = idmap.find(idb);
		if(miter == idmap.end())
		{
			miter = idmap.insert(pair<string, size_t>(idb, name_vec.size())).first;
			name_vec.push_back(idb);
			degree_vec.push_back(0);
		}
		j = miter->second;

		// lower triangle coding
		if(i>j)
			inx = i*(i-1)/2 + j;
		else
			inx = j*(j-1)/2 + i;

		if(included.find(inx) != included.end())
		{
			cerr << "Warning: Jump duplicated edge " << ida << ',' << idb << endl;
			continue;
		}

		if(compact_mode){
			fout.write((const char*)(&i), sizeof(USHRT));
			fout.write((const char*)(&j), sizeof(USHRT));
		}else{
			fout << ida << '\t' << idb << '\n';
		}

		degree_vec[i]++;
		degree_vec[j]++;

		included.insert(inx);
	}

	fin.close();
	fout.close();

	Nvertex = name_vec.size();
	Nedge = included.size();
	included.clear();

	if(compact_mode)
	{
		check(Nvertex > USHRT_MAX, "Error: Number of vertices exceeds " << USHRT_MAX << ".\nCannot use compact mode.")

		fout.open((output + ".map").c_str());
		for(i=0; i<Nvertex; i++) fout << name_vec[i] << '\n';
		fout.close();
	}

	adjacency_list = new size_t*[Nvertex];
	adjacency_list_edge_index = new size_t[2*Nedge];
	degree_cur = new size_t[Nvertex];
	candidate_arr = new size_t[Nvertex];
	flag = new bool[Nvertex];

	for (i=j=0 ; i<Nvertex ; i++)
	{
		candidate_arr[i] = i;

		// build up network structure
		adjacency_list[i] = adjacency_list_edge_index + j;
		j += degree_vec[i];
	}

	srand(0);

	for(i=0;i<Nrand;i++)
	{
		cout << "randomize " << i << endl;

		oss << output << '.' << i;

		max_ratio = ratio = 0;

		for(j=0;j<Nretry;j++)
		{
			random_shuffle(candidate_arr, candidate_arr + Nvertex);
			Nconnected = randomize(degree_vec, candidate_arr, adjacency_list, degree_cur, flag);

			ratio = (double)Nconnected/Nedge;

			if(ratio > max_ratio)
			{
				output_network(name_vec, degree_cur, adjacency_list, oss.str(), compact_mode);
				max_ratio = ratio;
			}

			if(max_ratio >= rthres) break;
		}

		if(max_ratio < rthres) cerr << "Warning: Insufficient connection ratio " << max_ratio << endl;
		oss.str(""); oss.clear();
	}

	delete[] adjacency_list;
	delete[] adjacency_list_edge_index;
	delete[] degree_cur;
	delete[] candidate_arr;
	delete[] flag;

	return 0;
}
