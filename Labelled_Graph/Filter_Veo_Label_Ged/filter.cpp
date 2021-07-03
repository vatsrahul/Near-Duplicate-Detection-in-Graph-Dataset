#include "veo.h"
#include "Gsim/GED.h"


using namespace std;

// For parsing the input graph dataset
void parseGraphDataset(ifstream &inp, vector<Graph> &graph_dataset, int &dataset_size);

// Sorts vertex and edge set of each graph in the dataset
void sortGraphDataset(vector<Graph> &graph_dataset);

// Graph comparator
bool graphComp(Graph &g1, Graph &g2);

// Returns the time from start to end in milliseconds
unsigned long long int clocksTosec(chrono::high_resolution_clock::time_point start, chrono::high_resolution_clock::time_point end);

// Displays the memory used by the program(in MB)
double memoryUsage();

// prints correct usage of program in case of an error
void usage();

// loose:   ./filter inp_file 1 simScore_threshold dataset-size res-file
// strict:  ./filter inp_file 2 simScore_threshold mismatch dataset-size res-file
//		     		  false/true : 0/1
// static:  ./filter inp_file 3 simScore_threshold mismatch noofbuckets dataset-size res-file
// dynamic: ./filter inp_file 4 simScore_threshold mismatch noofbuckets dataset-size res-file


void printingAndWritingInitialStatistics(int choice,double simScore_threshold,int dataset_size,const string res_dir,bool mismatch,int no_of_buckets)
{
	cout << "GSimJoin: VEO Similarity(filters)" << endl;
	cout << "Choice: " << choice << endl;
	cout << "Similarity Score Threshold: " << simScore_threshold << endl;
	cout << "Dataset Size: " << dataset_size << endl;

	ofstream stat_file(res_dir+"/stat_file.txt");
	stat_file << "GSimJoin: VEO Similarity(filters)" << endl;
	stat_file << "Choice: " << choice << endl;
	stat_file << "Similarity Score Threshold: " << simScore_threshold << endl;
	stat_file << "Dataset Size: " << dataset_size << endl;
	if(choice >= 2)
	{
		cout << "Mismatch: " << mismatch << endl;
		cout << "No of Buckets: " << no_of_buckets << endl;
		stat_file << "Mismatch: " << mismatch << endl;
		stat_file << "No of Buckets: " << no_of_buckets << endl;
	}
	stat_file.close();
}

void printingAndWritingFinalStatistics(double simScore_threshold,int choice,unsigned long looseCount,unsigned long strictCount,unsigned long VertexLabelFilterCount,bool isBucket,unsigned long EdgeLabelFilterCount,unsigned long SuffixFilterCount,bool mismatch,unsigned long mismatchCount,unsigned long simPairCount,int totalTimeTaken,const string res_dir,vector<long long int>& global_score_freq,unordered_map<unsigned, vector<pair<unsigned, double>>>& g_res)
{
    // Displaying stat file...
	unsigned CandidatePairs;
	if(choice >= 1)
	{
		cout << "Loose Filter Count: " << looseCount << endl;
		CandidatePairs = looseCount;
	}
	if(choice >= 2)
	{
		cout << "Strict Filter Count: " << strictCount << endl;
		CandidatePairs = strictCount;
	}
	if(choice >= 3)
	{
		cout << "Vertex Label Filter Count: " << VertexLabelFilterCount << endl;
		CandidatePairs = VertexLabelFilterCount;
	}
	if(choice >= 4)
	{
		cout << "Edge Label Filter Count: " << EdgeLabelFilterCount << endl;
		CandidatePairs = EdgeLabelFilterCount;
	}
	if(choice >= 5)
	{
		cout << "Suffix Filter Count: " << SuffixFilterCount << endl;
		CandidatePairs = SuffixFilterCount;
	}
	
	if(mismatch)
		cout << "Mismatch Filter Count: " << mismatchCount << endl;

	cout << "Final Similar Pair Count: " << simPairCount << endl;
	cout << "Memory used: " << memoryUsage() << " MB" << endl;
	cout <<"Total Time Taken: "<< totalTimeTaken << " milliseconds" << endl << endl	;
        
        ofstream stat_file(res_dir+"/stat_file.txt");
	stat_file.open(res_dir+"/stat_file.txt", ios::app);
	// Writing counts to stat file
	if(choice >= 1)
		stat_file << "Loose Filter Count: " << looseCount << endl;
	if(choice >= 2)
		stat_file << "Strict Filter Count: " << strictCount << endl;
	if(choice >= 3)
	{
		stat_file <<  "Vertex Label Filter Count: " << VertexLabelFilterCount << endl;
	}
	if(choice >= 4)
	{
		stat_file << "Edge Label Filter Count: " << EdgeLabelFilterCount << endl;
	}
	if(choice >= 5)
	{
		stat_file << "Suffix Filter Count: " << SuffixFilterCount << endl;
	}
	
	if(mismatch)
		stat_file << "Mismatch Filter Count: " << mismatchCount << endl;
	stat_file << "Final Similar Pair Count: " << simPairCount << endl;

	stat_file << "Memory used: " << memoryUsage() << " MB" << endl;
	stat_file <<"Total Time Taken: "<< totalTimeTaken << " milliseconds" << endl;
	stat_file.close();
	
	ofstream freq_file("./"+res_dir+"/freq_distr_file.txt");
	// for simScore==0
	freq_file << "0 " << global_score_freq[0] << endl; 
	for(int i=1; i<101; i++)
		freq_file << i << " " << global_score_freq[i] << endl;
	// for simScore==100
	freq_file << "101 " << global_score_freq[101] << endl; 
	freq_file.close();

	ofstream all_graph_file("./"+res_dir+"/all_graph_file.txt");	
	// Writing the result-set for each graph to the file for each graph
	for(auto g1 = g_res.begin(); g1 != g_res.end(); g1++)
	{
		for(auto g2 = g_res[g1->first].begin(); g2 != g_res[g1->first].end(); g2++)
		{
			all_graph_file << g1->first << " " << g2->first << " " << g2->second << endl;
		}
	}
	all_graph_file.close();



	// writing results in csv file

	ofstream stats("csv_file.csv", ios::app);
	//stat_file.open("csv_file.txt", ios::binary | std::ios_base::app);
	
	stats << totalTimeTaken << ", " <<CandidatePairs << ", " << ", ";
	if(simScore_threshold == 100)
		stats <<"\n\n";

	//stat_file.close();

}

int min_size=100000;  // min graph size in dataset.  Usd to convert ged to threshold

	unordered_set<unsigned> v_label;
	unordered_set<unsigned> e_label;
		// these 2 to count unique labels across all graph dataset

int main(int argc, char const *argv[])
{
	if(argc<6)
		usage();

	vector<Graph> graph_dataset; // input graph dataset

	// applying mismatch filter 
	bool mismatch=false; 
	// no. of buckets used in dynamic filter
	int no_of_buckets=0; 
	// true if no. of buckets is greater than 0 
	bool isBucket=false; 

	int choice = stoi(argv[2]);
	
	// Verifying args
	if(choice==1)
	{
		if(argc!=6)
			usage();
	}
	else if(choice==2)
	{
		if(argc!=7)
			usage();
		mismatch = (stoi(argv[argc-3])==1);
	}
	else if(choice > 2)
	{
		if(argc!=8)
			usage();
		mismatch = (stoi(argv[argc-4])==1);
		no_of_buckets = stoi(argv[argc-3]);
		isBucket = (no_of_buckets > 0);
	}
	else
		usage();

	int GED = stod(argv[3]); // similarity threshold 
	int dataset_size = stoi(argv[argc-2]); // size of input dataset
	const string res_dir = argv[argc-1]; // directory in which all stat files would be stored 
	mkdir(res_dir.c_str(),0777);

	ifstream dataset_file(argv[1]);
	if(!dataset_file.is_open())
	{
		cerr << "Unable to open dataset file" << endl;
		exit(0);
	}
	// parsing input graph-dataset
	parseGraphDataset(dataset_file, graph_dataset, dataset_size);
	cout << "Graph Dataset parsed.\n";

	// Sorts vertex and edge set of each graph in the dataset
	sortGraphDataset(graph_dataset);
	// sorts the graph dataset
	sort(graph_dataset.begin(), graph_dataset.end(), graphComp);
	cout << "Graph Dataset sorted.\n";

	unsigned long looseCount = 0; // No. of graphs filtered by loose size filter
	unsigned long strictCount = 0; // No. of graphs filtered by strict size filter
	unsigned long VertexLabelFilterCount = 0; // No. of graphs filtered by static index filter
	unsigned long EdgeLabelFilterCount = 0; // No. of graphs filtered by static index filter
	unsigned long PrefixFilterCount = 0; // No. of graphs filtered by static index filter
	unsigned long PositioningFilterCount = 0; // No. of graphs filtered by dynamic index filter
	unsigned long mismatchCount = 0; // No. of graphs filtered by mismatching filter
	unsigned long SuffixFilterCount = 0; // No. of graphs filtered by partition filter
	unsigned long simPairCount = 0; // No. of graph pairs having similarity score > threshold
	bool out = false; // a flag used to indicate whether graph is pruned or not 
	double simScore; // similarity score 

//////////////////////////////////////////////////////////////////////////////////
cout<<GED<<" "<<min_size<<" @\n";
	double simScore_threshold = (1.0*min_size)/(2*(1.0*min_size)+GED)*200.0;   ///////////// Converting GED to threshold
	simScore_threshold = floor(simScore_threshold);
//////////////////////////////////////////////////////////////////////////////////


	printingAndWritingInitialStatistics(choice,simScore_threshold,dataset_size,res_dir,mismatch,no_of_buckets);
	
	

	VEO veo_sim = VEO(simScore_threshold);
/*
	// static/dyanmic partition filter
	if(choice >= 3) 
	{
		veo_sim.index(graph_dataset, choice, isBucket, no_of_buckets); // index input graphs
		//veo_sim.calculate_sparse_table(graph_dataset, 0);
	}
	if(choice == 5)
		veo_sim.Preprocess_Suffix(graph_dataset, no_of_buckets);	// preprocessing for suffix filter
*/

	// Result-set for each graph as vector of other graph's gid and their similarity score as double
	unordered_map<unsigned, vector<pair<unsigned, double>>> g_res; // stores graph pair with similarity score 
	// Freq of simScore with range of 1% 0-1, 1-2, 1-3, ... 99-100% 
	vector<long long int> global_score_freq(102, 0); // stores sim-score frequency distribution of the dataset

 	// timestamping start time
	chrono::high_resolution_clock::time_point cl0 = chrono::high_resolution_clock::now();

////////////////////////////////////////////////////////////
	Preprocess_GED(argv[1], 2, GED);       /////// Preprocessing for GED with q-gram=4 and tau=5
////////////////////////////////////////////////////////////
	//unordered_map<unsigned, unordered_map<unsigned, unsigned> > temp;
	//temp.clear();


	for(int g1 = 0; g1 < graph_dataset.size(); g1++)
	{
		// size of current graph g1
		long double currSize = graph_dataset[g1].vertexCount + graph_dataset[g1].edgeCount; 
		//loose bound of PrevSize
		long double minPrevSize = floor(currSize/(long double)veo_sim.ubound);

	//	if(choice >= 3) 
	//		veo_sim.calculate_sparse_table(graph_dataset, g1, minPrevSize);	// this sparse table for g1 will be used in prefix and positioning filters.
	//	bool done = true;
		for(int g2 = g1-1; g2 >= 0; g2--)
		{
			double common = 0;
			out = false;
			// size of current graph g2
			long double PrevSize = graph_dataset[g2].vertexCount + graph_dataset[g2].edgeCount; 

			// loose filter
			if(PrevSize >= minPrevSize)	
				looseCount++;
			else
				break;

			if(choice >= 2)// strict filter
			{
				double maxIntersection = min(graph_dataset[g1].vertexCount, graph_dataset[g2].vertexCount) + min(graph_dataset[g1].edgeCount, graph_dataset[g2].edgeCount);
				// strict bound
				double strictBound = (double)200.0*maxIntersection/(currSize+PrevSize); 

				//strict filter
				if(simScore_threshold <= strictBound) 
					strictCount++;
				else
					continue;
			}
			if(choice >= 3){
				if(!out)
				{
					//if(done)
						//veo_sim.calculate_sparse_table(graph_dataset, g1),	// this sparse table for g1 will be used in prefix and positioning filters.
						//done = false;
					//out = veo_sim.PrefixFilter(graph_dataset[g1], graph_dataset[g2], g1, g2, choice, isBucket, no_of_buckets, PrefixFilterCount, simScore_threshold);
					out = veo_sim.VertexLabelFilter(graph_dataset[g1], graph_dataset[g2], g1, g2, v_label, e_label, VertexLabelFilterCount, simScore_threshold);
				}
			}
			if(choice >= 4){
				if(!out)
					//out = veo_sim.PositioningFilter(graph_dataset[g1], graph_dataset[g2], g1, g2, choice, isBucket, no_of_buckets, PositioningFilterCount, simScore_threshold);
					out = veo_sim.EdgeLabelFilter(graph_dataset[g1], graph_dataset[g2], g1, g2, v_label, e_label, EdgeLabelFilterCount, simScore_threshold);
				
			}
			
			if(choice == 5) // suffix filter
			{ 

				if(!out)
				out = veo_sim.SuffixFilter(graph_dataset[g1], graph_dataset[g2], g1, g2, simScore_threshold, isBucket, no_of_buckets, SuffixFilterCount);
			}
			
			/*if(!out){ // vertex filter

				out = veo_sim.VertexFilter(graph_dataset[g1], graph_dataset[g2], g1, g2, simScore_threshold);
				if(!out)
					suffixCount++;
			}*/

			if(out)
				continue;
			else if(mismatch) 
			{
				// mismatching filter
				out = veo_sim.mismatchingFilter(graph_dataset[g1], graph_dataset[g2], common, simScore_threshold);
				if(!out)
					mismatchCount++;
			}
			if(!out)
			{

		

				// naive computation of VEO similarity
				simScore = veo_sim.computeSimilarity(graph_dataset[g1], graph_dataset[g2], common, v_label, e_label);
				if(simScore >= simScore_threshold)
				{//cout<<graph_dataset[g1].gid<<" "<<graph_dataset[g2].gid<<" ";
				//simPairCount++;
					int edit_dist = Get_GED(graph_dataset[g2].gid, graph_dataset[g1].gid);//cout<<edit_dist<<"\n";
				if(edit_dist > GED)
					continue;
					//temp[graph_dataset[g2].gid][graph_dataset[g1].gid]=2;
					//temp[graph_dataset[g1].gid][graph_dataset[g2].gid]=2;
					
					// Incrementing count... 
					if(simScore==0.0)
					{
						// Disjoint graphs
						global_score_freq[0]++;
					}
					else if(simScore==100.0)
					{
						// Identical graphs
						global_score_freq[101]++;
					}
					else
					{
						// example: 54.5% will be mapped to index 60
						// example: 0.5% will be mapped to index 1
						// example: 99.5% will be mapped to index 100
						global_score_freq[(int)ceil(simScore)]++;
					}
					g_res[graph_dataset[g1].gid].push_back(make_pair(graph_dataset[g2].gid, simScore));
					simPairCount++;
				}
			}
		}
	}
	/* ifstream ifp("Gsim/2k-gsim45.txt");
  
int i,j,ged,cot=0;
 long double mx=0.0, mn = 1.0,sum=0.0;

          
    while(ifp>>i and  ifp>>j and ifp>>ged){
		//if((temp.find(i)!=temp.end() && temp[i],find(j)!=temp[i].end()))
		if(temp[i][j]==2 and temp[j][i]==2)
		{ cot++;}
		else{
			cout<<i<<" "<<j<<"\n";
			cout<<"not found\n";
			
		}
            
        }cout<<"all matche "<<cot<<endl;*/

 	////////////////////////////////////////////////////////
	 Clean_GED();							///////////////Cleaning the mem alloted.
	///////////////////////////////////////////////////////

	 // timestamping end time
	chrono::high_resolution_clock::time_point cl2 = chrono::high_resolution_clock::now();
	int totalTimeTaken = (clocksTosec(cl0,cl2));

    printingAndWritingFinalStatistics(simScore_threshold,choice,looseCount,strictCount,VertexLabelFilterCount,isBucket,EdgeLabelFilterCount,SuffixFilterCount,mismatch,mismatchCount,simPairCount,totalTimeTaken,res_dir,global_score_freq,g_res);

	return 0;
}

// parses the input graph dataset and query graph
void parseGraphDataset(ifstream &inp, vector<Graph> &graph_dataset, int &dataset_size)
{
	int size;
	inp >> size;
	if(dataset_size == -1)
		dataset_size=size;
	graph_dataset.resize(dataset_size);

	for(auto g_iter = graph_dataset.begin(); g_iter != graph_dataset.end(); g_iter++)
	{
		g_iter->readGraph(inp, v_label, e_label);
		min_size = min(min_size, (int)g_iter->vertexCount+(int)g_iter->edgeCount);
	}
}

bool graphComp(Graph &g1, Graph &g2)
{
	return g1.vertexCount+g1.edgeCount < g2.vertexCount+g2.edgeCount;
}

// Sorts vertex and edge set of graph dataset
void sortGraphDataset(vector<Graph> &graph_dataset) 
{
	for(int i = 0; i < graph_dataset.size(); i++)
	{
		sort(graph_dataset[i].vertices.begin(), graph_dataset[i].vertices.end()); // sort vertex-set 
		sort(graph_dataset[i].edges.begin(), graph_dataset[i].edges.end()); // sort edge-set
	}
}

// Returns the time from start to end in seconds
unsigned long long int clocksTosec(chrono::high_resolution_clock::time_point start, chrono::high_resolution_clock::time_point end)
{
	return (unsigned long long int)(1e-6*chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

// Displays the memory used by the program(in MB)
double memoryUsage()
{
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	return r_usage.ru_maxrss/1024.0;
}

// prints correct usage of program in case of an error
void usage(){
	cerr << "usage: ./filter input_file choice simScore-threshold dataset-size res-dir" <<endl;

	cerr << "Available choices: " << endl;
	cerr << endl;

	cerr << "1 loose : 1" << endl;
	cerr << endl;

	cerr << "2 loose + strict 			  : 2, mismatch=false noofbuckets=0" << endl;
	cerr << "3 loose + strict + mismatch  : 2, mismatch=true  noofbuckets=0" << endl;
	cerr << endl;

	cerr << "4 loose + strict + static 				: 3, mismatch=false noofbuckets = 0" << endl;
	cerr << "5 loose + strict + static + 2 buckets  : 3, mismatch=false noofbuckets = 2" << endl;
	cerr << "6 loose + strict + static + 5 buckets  : 3, mismatch=false noofbuckets = 5" << endl;
	cerr << "7 loose + strict + static + 10 buckets : 3, mismatch=false noofbuckets = 10" << endl;
	cerr << "8 loose + strict + static + 10 buckets : 3, mismatch=true  noofbuckets = 10" << endl;
	cerr << "9  loose + strict + dynamic 			  :  4, mismatch=false noofbuckets = 0" << endl;
	cerr << "10 loose + strict + dynamic + 2 buckets  :  4, mismatch=false noofbuckets = 2" << endl;
	cerr << "11 loose + strict + dynamic + 5 buckets  :  4, mismatch=false noofbuckets = 5" << endl;
	cerr << "12 loose + strict + dynamic + 10 buckets :  4, mismatch=false noofbuckets = 10" << endl;
	cerr << "13 loose + strict + dynamic + 10 buckets :  4, mismatch=true  noofbuckets = 10" << endl;
	cerr << endl;

	cerr << "loose:   ./filter inp_file 1 simScore_threshold dataset-size res-dir\n";
	cerr << "strict:  ./filter inp_file 2 simScore_threshold mismatch dataset-size res-dir\n";
//		     											  false/true : 0/1
	cerr << "static:  ./filter inp_file 3 simScore_threshold mismatch noofbuckets dataset-size res-dir\n";
	cerr << "dynamic: ./filter inp_file 4 simScore_threshold mismatch noofbuckets dataset-size res-dir\n";
	exit(0);
}

/*
#include "veo.h"

using namespace std;

// For parsing the input graph dataset
void parseGraphDataset(ifstream &inp, vector<Graph> &graph_dataset, int &dataset_size);

// Sorts vertex and edge set of graph dataset
void sortGraphDataset(vector<Graph> &graph_dataset);

// Returns the time from start to end in seconds
unsigned long long int clocksTosec(chrono::high_resolution_clock::time_point start, chrono::high_resolution_clock::time_point end);

// Displays the memory used by the program(in MB)
double memoryUsage();

// prints correct usage of program in case of an error
void usage(); 

// $ ./naive inp-file simScore_threshold dataset-size res-file

int min_size=100000;  // min graph size in dataset.  Usd to convert ged to threshold
unordered_set<unsigned> v_label;
unordered_set<unsigned> e_label;
// these 2 to count unique labels across all graph dataset

int main(int argc, char const *argv[])
{
	if(argc!=5)
		usage();

	int GED = stod(argv[2]);  // threshold to write only those graph pairs to all_graph_file.txt
	int dataset_size = stoi(argv[3]); // size of input dataset
	const string res_dir = argv[4]; // directory in which all stat files would be stored
	//cout << dataset_size;
	mkdir(res_dir.c_str(),0777);
	mkdir((res_dir+"/graph_details/").c_str(),0777);
	vector<Graph> graph_dataset; // input graph dataset

	ifstream dataset_file(argv[1]);
	if(!dataset_file.is_open())
	{
		cerr << "Unable to open dataset file" << endl;
		exit(0);
	}

	// parsing input dataset 
	parseGraphDataset(dataset_file, graph_dataset, dataset_size);
	cout << argv[1] << ": Graph Dataset parsed." << endl;

	//cout << 1;


	sortGraphDataset(graph_dataset); // to sort vertex and edge set
	cout << "All graphs in dataset sorted." << endl;
    
	

	ofstream all_graph_file("./"+res_dir+"/all_graph_file.txt");
	all_graph_file.close();

    double simScore_threshold = (1.0*min_size)/(2*(1.0*min_size)+GED)*200.0;   ///////////// Converting GED to threshold
	simScore_threshold = floor(simScore_threshold);

	// Result-set for each graph as vector of other graph's gid and their similarity score as double
	vector<pair<unsigned, double>> g_res; // stores graph pair with the score of a specific graph
	vector<unsigned long long int> g_time(graph_dataset.size()); // stores time required for each graph in the dataset
	unsigned long long int global_time = 0; // total time taken for similarity computation 
	// Freq of simScore with range of 1% 0-1, 1-2, 1-3, ... 99-100% 
	vector<int> score_freq(102,0); // stores sim-score frequency distribution of a particular graph
	vector<long long int> global_score_freq(102, 0); // stores sim-score frequency distribution of the dataset

	double simScore; // similarity score
	unsigned long long int simPairCount = 0; // no. of graph pairs having similarity score > threshold

 	// timestamping start time
	chrono::high_resolution_clock::time_point cl0 = chrono::high_resolution_clock::now();

	// For clock-time calculation
	chrono::high_resolution_clock::time_point clTemp0, clTemp1; 


//unordered_map<unsigned, unordered_map<unsigned, unsigned> > temp;
	//temp.clear();

	for(int g1 = 1; g1<graph_dataset.size(); g1++)
	{
		clTemp0 = chrono::high_resolution_clock::now();

		for(int g2 = g1-1; g2 >= 0; g2--)
		{

			// Similarity Calculation...
			double common = 0;

			//cout << 1;
			double simScore = computeSimilarity(graph_dataset[g1], graph_dataset[g2], common, v_label, e_label);
			
			// Incrementing count... 
			if(simScore==0.0)
			{
				// Disjoint graphs
				score_freq[0]++;
				global_score_freq[0]++;
			}
			else if(simScore==100.0)
			{
				// Identical graphs
				score_freq[101]++;
				global_score_freq[101]++;
			}
			else
			{
				// example: 54.5% will be mapped to index 60
				// example: 0.5% will be mapped to index 1
				// example: 99.5% will be mapped to index 100
				score_freq[(int)ceil(simScore)]++;
				global_score_freq[(int)ceil(simScore)]++;
			}

			// Storing only those graph pairs which have similarity above (simScore_threshold)%
			if(simScore>=simScore_threshold)
			{
				//temp[graph_dataset[g2].gid][graph_dataset[g1].gid]=2;
				//temp[graph_dataset[g1].gid][graph_dataset[g2].gid]=2;
			
				g_res.push_back(make_pair(graph_dataset[g2].gid, simScore));
				simPairCount++;
			}
		}

		clTemp1 = chrono::high_resolution_clock::now();
		g_time[g1] = clocksTosec(clTemp0,clTemp1); // graph's similarity calculation time
		global_time += g_time[g1]; // dataset's similarity calculation time

		// Creating Result Files for graph g1
		ofstream gfile;
		all_graph_file.open("./"+res_dir+"/all_graph_file.txt",ios::app);
		gfile.open("./"+res_dir+"/graph_details/g_"+to_string(g1)+"_"+to_string(graph_dataset[g1].gid)+"_sim.txt");

		// Writing the result-set for each graph to the file for each graph
		gfile << g_res.size() << endl;
		for(auto g_iter = g_res.begin(); g_iter != g_res.end(); g_iter++)
		{
			gfile << graph_dataset[g1].gid << " " << g_iter->first << " " << g_iter->second << endl;
			all_graph_file << graph_dataset[g1].gid << " " << g_iter->first << " " << g_iter->second << endl;
		}
		g_res.clear();
		all_graph_file.close();

		// Writing g1's simScore-freq
		// for simScore==0
		gfile << "0 " << score_freq[0] << endl;
		score_freq[0] = 0;
		for(int i=1; i<101; i++)
		{
			gfile << i << " " << score_freq[i] << endl;
			score_freq[i] = 0;
		}
		// for simScore==100
		gfile << "101 " << score_freq[101] << endl; 
		score_freq[101] = 0;
		gfile << g_time[g1] << endl;
		gfile << global_time << endl;
		gfile.close();
	}
 	
	 ifstream ifp("2k_gsim23.txt");  
	unsigned i,j,ged,cot=0;
    while(ifp>>i and  ifp>>j and ifp>>ged){
		//if((temp.find(i)!=temp.end() && temp[i],find(j)!=temp[i].end()))
		if(temp[i][j]==2 and temp[j][i]==2)
		{ cot++;}
		else{
			cout<<i<<" "<<j<<"\n";
			cout<<"not found\n";
		}
            
        }cout<<"all matche "<<cot<<endl;

	 
	 
	 // timestamping end time
	chrono::high_resolution_clock::time_point cl1=chrono::high_resolution_clock::now();	
	
	cout << "GSimJoin: VEO Similarity(naive)" << endl;
	cout << "Dataset size: " << dataset_size << endl;
	cout << "Similarity Score Threshold: " << simScore_threshold << endl;
	cout << "Similar Graph Pairs: " << simPairCount << endl;
	cout << "Memory used: " << memoryUsage() << " MB" << endl;
	cout << "Similarity Time: "<< global_time << " milliseconds" << endl;
	cout << "Total Time Taken: "<< (clocksTosec(cl0,cl1))  << " milliseconds" << endl;

	ofstream stat_file("./"+res_dir+"/stat_final.txt");
	stat_file << "GSimJoin: VEO Similarity(naive)" << endl;
	stat_file << "Dataset size: " << dataset_size << endl;
	stat_file << "Similarity Score Threshold: " << simScore_threshold << endl;
	stat_file << "Similar Graph Pairs: " << simPairCount << endl;
	stat_file << "Memory used: " << memoryUsage() << " MB" << endl;
	stat_file << "Similarity Time: "<< global_time << " milliseconds" << endl;
	stat_file << "Total Time Taken: "<< (clocksTosec(cl0,cl1))  << " milliseconds" << endl;
	stat_file.close();
	
	ofstream freq_file("./"+res_dir+"/freq_distr_file.txt");
	// for simScore==0
	freq_file << "0 " << global_score_freq[0] << endl; 
	for(int i=1; i<101; i++)
		freq_file << i << " " << global_score_freq[i] << endl;
	// for simScore==100
	freq_file << "101 " << global_score_freq[101] << endl; 
	freq_file.close();

	return 0;
}

// For parsing the input graph dataset
void parseGraphDataset(ifstream &dataset_file, vector<Graph> &graph_dataset, int &dataset_size)
{
	int size;
	dataset_file >> size;
	if(dataset_size == -1)
		dataset_size = size;
	graph_dataset.resize(dataset_size);
	for(auto g_iter = graph_dataset.begin(); g_iter != graph_dataset.end(); g_iter++)
	{	g_iter->readGraph(dataset_file, v_label, e_label);
		min_size = min(min_size, (int)g_iter->vertexCount+(int)g_iter->edgeCount);
	}	
}

// Sorts vertex and edge set of graph dataset

void sortGraphDataset(vector<Graph> &graph_dataset)
{
	for(int i = 0; i < graph_dataset.size(); i++)
	{
		//// sort vertex-set ON BASIS OF LABELS
		sort(graph_dataset[i].vertices.begin(), graph_dataset[i].vertices.end(),    
		[&] (const int a, const int b) {return graph_dataset[i].vid_to_vc[a] < graph_dataset[i].vid_to_vc[b];}); 
		
		// edges are sorted on bases of labels only as it contains nothing but label ascii values.
		sort(graph_dataset[i].edges.begin(), graph_dataset[i].edges.end()); // sort edge-set
	}
}

// Returns the time from start to end in seconds
unsigned long long int clocksTosec(chrono::high_resolution_clock::time_point start, chrono::high_resolution_clock::time_point end)
{
	return (unsigned long long int)(1e-6*chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

// Displays the memory used by the program(in MB)
double memoryUsage()
{
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	return r_usage.ru_maxrss/1024.0;
}

void usage()
{
	cerr << "usage: ./naive inp-file simScore_threshold dataset-size res-file" <<endl;
	exit(0);
}
*/