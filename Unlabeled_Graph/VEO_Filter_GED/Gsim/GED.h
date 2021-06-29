#include "common1.h"
#include "path_join1.h"
#include "path_search1.h"

// setting variables
char vf_order = '1', version = '0';
bool count_filter = true, min_edit_filter = true, local_lab_filter = true, deg_match_cond = true;
bool filter_only = false, print_ans = false, print_more = false;

// trivial data structure
int tau = 5, over_tau = 6, qs = 5, path_qs = 9, under_qs = 4;
unsigned max_vnum;
bool* done;

// main data structure
vector<Graph1> gdb;
vector<Path_gram*>* pdb;		// path grams, sorted in frequency
Hash_gram** hdb;				// hash values only, sorted in value
vector<unsigned> uid;			// underflow: ids of un-indexed graphs
set<unsigned>** eds;			// edit effect: ids of paths controlled by each vertex
unsigned* glen;                 // prefix length
fmap path_freq;
gmap path_graph;

// assistant data structure
unsigned gdb_size;				// number of data graphs
unsigned qry_size;				// number of query graphs
int* status;
unsigned global_esize;
unsigned* global_elabs;
unsigned** card;			    // size of each edit set
unsigned** tau_card;			// size of top-tau set
vector<int>* elab;				// edge labels
vector<unsigned>* ecnt;			// e-label counts
vector<int>* vlab;				// vertex labels
vector<unsigned>* vcnt;			// v-label counts
vector<iter_graph>* git;
Freq_comp fcomp(&path_freq);

unordered_map<unsigned,unsigned> name_to_index;

void Preprocess_GED(string p, int under_qs, int tau);
int Get_GED(int i, int j);
void Clean_GED();

void usage1() {
	cerr << "GSim - Graph similarity queries with path q-grams\n\n"
		<< "Usage: GSim [OPTION...] [DB]... [Q] [TAU]\n\n"
		<< "Examples:\n"
		<< "GSimJoin -m0 -o1 -t2 db        # Self-join on db with 4-grams & threshold 2, applying count, min-edit, local-label filters and optimized verification\n"
		<< "GSimJoin -m2 -o0 -t3 db q 3 1  # Similarity search on db using q with 3-grams & threshold 1, applying all filters\n" << endl
		<< "Options:\n"
		<< "  -f            execute filtering only without verification\n"
		<< "  -h            display this help information\n"
		<< "  -m[i]         specify work mode\n"
		<< "                [0 - self-join (default) | 1 - R-S join | 2 - search]\n"
		<< "  -o[i]         specify verification method\n"
		<< "                [0 - basic | 1 - optimized, available with -t2]\n"
		<< "  -t[i]         apply filter combination\n"
		<< "                [0 - count filter | 1 - + min-edit filter | 2 - + local-label filter | 3 - + degree-based match (default)]\n"
		<< "  -v            verbosely print results with cerr (may sacrifice runtime performance)\n" << endl
		<< "Arguments:\n"
		<< "  DB            data graph file (mandatory)\n"
		<< "  Q             q-gram size - positive number of hops (default 4)\n"
		<< "  TAU           non-negative GED threshold (default 2)\n\n"
		<< "*This* GSim defaults to:\n"
		<< "-m0 -o1 -t3 [gdb] 4 2\n\n"
		<< "Report bugs to <xiangzhao@nudt.edu.cn>." << endl;
	exit(0);
}

void parse_data(char* file) {
	ifstream fin(file, ifstream::in);
	if (!fin) {
		cerr << file << " is not open or found!" << endl;
		exit(1);
	} else {
		cout << "Loading...\n" << endl;
	}

	done = init_bool_array(true, LAG_SIZE);

	set<int> edge_label;
	char tag;
	int temp;
	fin >> temp;  // ignored value. dataset sizze;
	fin >> tag;
	while (!fin.eof()) {
		gdb.resize(gdb.size() + 1);
		tag = gdb.back().read(fin, tag, edge_label);
		if (gdb.back().size() <= tau || !gdb.back().connectivity() 
			|| !gdb.back().simplicity()) {	// discard disconnect graphs
			cout<<gdb.back().name<<" ";
			if(gdb.back().size()<=tau){
				cout<<"a";
			}
			cout<<"\n";
			gdb.resize(gdb.size() - 1);
		}
	}
	fin.close();

	global_esize = edge_label.size();
	global_elabs = new unsigned[global_esize];
	unsigned i = 0;
	for (set<int>::const_iterator it = edge_label.begin(), end = edge_label.end(); it != end; ++ i, ++ it) {
		global_elabs[i] = *it;
	}

	sort(gdb.begin(), gdb.end(), Graph_comp());	// sort graphs in ascending order w.r.t. |G|=|v|+|E|
	gdb_size = gdb.size();
	max_vnum = gdb.back().vertex_num;

int c=0;
	for(auto g: gdb){
		name_to_index[g.name]=c++;
	}

#ifdef UNIX
	path_freq.set_empty_key(0);
	path_graph.set_empty_key(0);
	path_graph.set_deleted_key(OVER_BOUND);
#endif
}

void load_query(char* file) {
	ifstream fin(file, ifstream::in);
	if (!fin) {
		cerr << file << " is not open or found!" << endl; exit(1);
	}

	char tag;
	fin >> tag;
	while (!fin.eof()) {
		gdb.resize(gdb.size() + 1);
		tag = gdb.back().read(fin, tag);
		/*if (!gdb.back().connectivity()) {
			gdb.resize(gdb.size()-1);
		} else {*/
			if (gdb.back().vertex_num > max_vnum) max_vnum = gdb.back().vertex_num; 
		//}
	}
	fin.close();
}

void init_mem()
{
	unsigned total_size = gdb.size();
	qry_size = total_size - gdb_size;
	pdb = new vector<Path_gram*>[total_size];
	hdb = new Hash_gram*[total_size];
	glen = new unsigned[total_size];
	memset(glen, 0, sizeof(unsigned)*total_size);
	eds = new set<unsigned>*[total_size];
	card = new unsigned*[total_size];
	tau_card = new unsigned*[total_size];
	for (unsigned i = 0; i != total_size; ++ i) {
		eds[i] = new set<unsigned>[gdb[i].vertex_num];
		card[i] = new unsigned[gdb[i].vertex_num];
		tau_card[i] = new unsigned[over_tau];
		memset(tau_card[i], 0, sizeof(unsigned)*over_tau);
		gdb[i].vlabel_list = new int[gdb[i].vertex_num];
		gdb[i].vid_to_vlabel = new int[gdb[i].vertex_num];
		gdb[i].eid_to_elabel = new vector<int>[gdb[i].vertex_num];
		gdb[i].elabel_list = new int[gdb[i].edge_num];
		for (unsigned j = 0; j != gdb[i].vertex_num; ++ j) {
			gdb[i][j].elabs = new unsigned[global_esize];
			memset(gdb[i][j].elabs, 0, global_esize*sizeof(unsigned));
			gdb[i].eid_to_elabel[j].resize(gdb[i][j].degree);
		}
	}
	status = new int[total_size];
	memset(status, 0, sizeof(int)*total_size);	// default 0 (class-1), can be 1 (class-2), -1 (class-3)
	elab = new vector<int>[total_size];
	ecnt = new vector<unsigned>[total_size];
	vlab = new vector<int>[total_size];
	vcnt = new vector<unsigned>[total_size];
	git = new vector<iter_graph>[gdb_size];	// only used in join
}

void Clean_GED() {
	cout << "\nDone..." << endl;
	delete[] done;
	for (unsigned i = 0; i != gdb.size(); ++ i) {
		/*for (unsigned j = 0; j != pdb[i].size(); ++ j) {
			delete pdb[i][j];
		}
		for (unsigned j = 0; j != gdb[i].vertex_num; ++ j) {
			delete[] gdb[i][j].elabs;
		}
		delete[] eds[i];
		delete[] hdb[i];
		delete[] card[i];
		delete[] tau_card[i];*/
	}
	delete[] global_elabs;
	delete[] pdb;
	delete[] eds;
	delete[] hdb;
	delete[] glen;
	delete[] status;
	delete[] card;
	delete[] tau_card;
	delete[] elab;
	delete[] ecnt;
	delete[] vlab;
	delete[] vcnt;
	delete[] git;
}

void Preprocess_GED(string p, int under_qs_local=4, int tau_local=5){
	char* pchar;
	bool file_set = false, q_set = false, tau_set = false;
	char dfile[256] = "gdb", qfile[256], comb;

				strcpy(dfile, p.c_str());

				under_qs = under_qs_local;
				qs = under_qs + 1;
				path_qs = (qs << 1) - 1;
				q_set = true;

				tau = tau_local;
				over_tau = tau + 1;
				tau_set = true;


				parse_data(dfile);
	
	cout << "Setting:";
	switch (version) {
		case '0': 
			cout << " Self-join -"; break;
		case '1': 
			if (strlen(qfile) == 0) { cerr << "The other data file not provided!\n" << endl; usage1(); }
			cout << " R-S join -"; break;
		case '2': 
			if (strlen(qfile) == 0) { cerr << "Query file not set!\n" << endl; usage1(); }
			cout << " Full-search -"; break;
		default:
			cerr << "No such outlet!\n" << endl;
			usage1();
	}
	if (count_filter) cout << " count ";
	if (min_edit_filter) cout << "min-edit ";
	if (local_lab_filter) cout << "label ";
	if (deg_match_cond) cout << "degree ";
	cout << "ft.";  // filters
	
	if (filter_only) cout << " ONLY";
	else {
		switch (vf_order) {
			case '0': cout << " + naive vf."; break;
			case '1': cout << " + improve vf."; break;  // verification
		}
	}


	init_mem();
			cout << "\n         Num = " << gdb_size << " Q = " << under_qs << " TAU = " << tau << " DB = " << dfile << endl;
			print_more = false;
			print_ans = true; 

//run_min_prefix();
//return;
	vectorize_label();
/*	generate_all_path();
	count_all_path();	
	join_min_prefix_index();

	iter_graph temp;
	for (iter_graph it = path_graph.begin(), end = path_graph.end(); it != end;) {
		if (it->second.size() < 2) {
			temp = it;
			++ it;
			path_graph.erase(temp);
		} else {
			for (set<unsigned>::const_iterator is = it->second.begin(), end = it->second.end(); is != end; ++ is) {
				git[*is].push_back(it);
			}
			++ it;
		}
	}*/
	//opt_order_join();

}

int Get_GED(int i, int j){

			int ans;
			ans=0;
			Priority* pri; 
			pri = new Priority(name_to_index[i],name_to_index[j]);
			//cout << pri->lgid << "(" << pri->lg->name << ") - " << pri->rgid << "(" << pri->rg->name << ") " ;
			
			int Edit_dist = compute_rud_dist(pri, ans, false);
			delete pri;
	return Edit_dist;
}
