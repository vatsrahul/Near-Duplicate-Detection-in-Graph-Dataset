#include "graph.h"
#include<unordered_set>
#include <map>

class VEO{
public:
	double ubound;
	map< pair<unsigned, unsigned>, unsigned long > rank;
	vector<vector<unsigned long>> rankList;
	vector<vector<pair<unsigned long ,unsigned long>>> InvertedIndex;
	unordered_map<int,pair<int,int>> sparse_table;    // a hash_table
	vector<vector<vector<unsigned>>> bucket;
	vector<vector<unsigned>> graph_bucket;
	VEO(double threshold)
	{
		ubound = double((double)(200.0/threshold) -1.0);
		cout << "ubound: " << ubound << endl;
	}
	bool PrefixFilter(Graph &g2, Graph &g1, int index2, int index1, int mode, bool isBucket, int no_of_buckets, long unsigned &indexCount, double threshold);
	bool PositioningFilter(Graph &g2, Graph &g1, int index2, int index1, int mode, bool isBucket, int no_of_buckets, long unsigned &indexCount, double threshold);
	void Preprocess_Suffix(vector<Graph> &graph_dataset, int no_of_buckets);
	bool SuffixFilter(Graph &g2, Graph &g1, int index2, int index1, double threshold, bool isBucket, int no_of_buckets, long unsigned &SuffixFilterCount);
	bool VertexFilter(Graph &g2, Graph &g1, int index2, int index1, double threshold);
	void ranking(vector<Graph> &graphDataset);
	void buildPrefix(vector<Graph> &graphDataset, int type, bool isBucket, int b);
	void calculate_sparse_table(vector<Graph> &graphDataset, int g_ind);
	bool mismatchingFilter(Graph &g1, Graph &g2, double &c, double threshold);
	void sortGraphDataset(vector<Graph> &graphDataset);
	double computeSimilarity(Graph &g1, Graph &g2, double &commonV);
	void sortVertexEdge(vector<Graph> &graphDataset);
	void printlist(int i);
	void index(vector<Graph> &graphDataset, int type, bool isBucket, int b);
};
