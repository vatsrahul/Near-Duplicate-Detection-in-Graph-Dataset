#include "veo.h"
#include <cmath>
#include <bits/stdc++.h>

bool freqComp(pair<pair<unsigned, unsigned>,unsigned> v1, pair<pair<unsigned, unsigned>,unsigned> v2)
{ 
    return (v1.second > v2.second); 
}

// Edge comparator
bool edgeComp(pair<unsigned, unsigned> &a, pair<unsigned, unsigned> &b,Graph &g1,Graph &g2)// since we have stored labels in edges, so we need not use vid_to_vc here.
{
	if(a.first == b.first)
		return a.second < b.second;
	return a.first < b.first;
}

// Intersection of vertex lists of 2 graphs
double intersection_vertices(vector<unsigned> &s1, vector<unsigned> &s2,unsigned gid1,unsigned gid2,Graph &g1,Graph &g2, unordered_set<unsigned>& v_label)
{
	double common = 0;

	for(auto label : v_label){
		common += min(g1.VertexLabelMap[label], g2.VertexLabelMap[label]);
	}


	/*unsigned s1_iter = 0;
	unsigned s2_iter = 0;
	
	while(s1_iter < s1.size() && s2_iter < s2.size())
	{
		if(g1.vid_to_vc[s1[s1_iter]] < g2.vid_to_vc[s2[s2_iter]])
			s1_iter++;
		else if(g1.vid_to_vc[s1[s1_iter]] > g2.vid_to_vc[s2[s2_iter]])
			s2_iter++;
		else
		{
			
			s1_iter++;
			s2_iter++;
			common++;		
		}
	}*/

	return common;
}

// Intersection of edge lists of 2 graphs
double intersection_edges(vector<unsigned> &s1, vector<unsigned> &s2,unsigned gid1,unsigned gid2,Graph &g1,Graph &g2, unordered_set<unsigned>& e_label)
{

	double common = 0;

	for(auto label : e_label){

		common += min(g1.EdgeLabelMap[label], g2.EdgeLabelMap[label]);
	}

/*	unsigned s1_iter = 0;
	unsigned s2_iter = 0;

	while(s1_iter < s1.size() && s2_iter < s2.size())
	{
		if(s1[s1_iter] < s2[s2_iter])
			s1_iter++;
		else if(s1[s1_iter] > s2[s2_iter])
			s2_iter++;
		else
		{
			
			s1_iter++;
			s2_iter++;
			common++;		
		}
	}
*/
	return common;
}

// VEO Similarity computation
double VEO:: computeSimilarity(Graph &g1, Graph &g2, double &commonV, unordered_set<unsigned>& v_label, unordered_set<unsigned>& e_label)
{
	if(commonV == 0)
		commonV = intersection_vertices(g1.vertices, g2.vertices,g1.gid,g2.gid,g1,g2, v_label);
	commonV+=intersection_edges(g1.edges, g2.edges,g1.gid,g2.gid,g1,g2, e_label);

	//cout << commonV << " "<<g1.vertexCount+g2.vertexCount+g1.edgeCount+g2.edgeCount<<"\n";

	double simScore = (double)(200.0*(commonV)/(double)(g1.vertexCount+g2.vertexCount+g1.edgeCount+g2.edgeCount));

	return simScore;
}


/////////////////////////////////////////////////////////////////////////////////////////////



// Sorts vertex and edge set of graph dataset
void VEO:: sortVertexEdge(vector<Graph> &graph_dataset)
{
	for(int i = 0; i < graph_dataset.size(); i++)
	{
		//sort(graph_dataset[i].vertices.begin(), graph_dataset[i].vertices.end());
		//sort(graph_dataset[i].edges.begin(), graph_dataset[i].edges.end());
	}
}

void VEO:: printlist(int k)
{
	for(int i = 0; i < rankList[k].size(); i++)
		cout << i  << ": "<< rankList[k][i] << endl;
}


unsigned long r; 	// to store max value of a rank/attribute ; to be used in suffix filter

//globaly ranks vertices and edges together
void VEO:: ranking(vector<Graph> &graph_dataset) 
{ return;
	// traversing the graph-dataset
/*	for(int g_ind = 0; g_ind < graph_dataset.size(); g_ind++)
	{
		// traversing the vertex-set

		for(int vtx_ind = 0; vtx_ind < graph_dataset[g_ind].vertices.size(); vtx_ind++)
		{
			pair<unsigned, unsigned> vtx_pair = make_pair(graph_dataset[g_ind].vertices[vtx_ind], graph_dataset[g_ind].vertices[vtx_ind]);
			if(rank.find(vtx_pair) != rank.end())
				rank[vtx_pair]++;
			else
				rank[vtx_pair] = 0;
		}
		// traversing the edge-set
		for(int edge_ind = 0; edge_ind < graph_dataset[g_ind].edges.size(); edge_ind++)
		{
			pair<unsigned, unsigned> edge_pair;
			if(graph_dataset[g_ind].edges[edge_ind].first > graph_dataset[g_ind].edges[edge_ind].second)
				edge_pair = make_pair(graph_dataset[g_ind].edges[edge_ind].second, graph_dataset[g_ind].edges[edge_ind].first);
			else
				edge_pair = make_pair(graph_dataset[g_ind].edges[edge_ind].first, graph_dataset[g_ind].edges[edge_ind].second);
			if(rank.find(edge_pair) != rank.end())
				rank[edge_pair]++;
			else
				rank[edge_pair] = 0;
		}
	}

	vector<pair<pair<unsigned, unsigned>, unsigned long> > freqList;
	copy(rank.begin(), rank.end(), back_inserter(freqList));
	sort(freqList.begin(), freqList.end(), freqComp);

	// ranking each vertex and each edge
	r = 1;
	for(auto entry = freqList.begin(); entry!=freqList.end(); entry++)
	{
		rank[entry->first] = r;
		r++;
	}
//	cout<<"Max rank is "<<r<<"\n";
	// setting the size of invertecindex list as max value of rank - r
	InvertedIndex.resize(r);
	*/

}

void VEO:: buildPrefix(vector<Graph> &graph_dataset, int mode, bool isBucket, int no_of_buckets)
{return ;/*
	rankList.resize(graph_dataset.size());
	if(isBucket)
		bucket.resize(graph_dataset.size());
	
	for(int g_ind = 0; g_ind < graph_dataset.size(); g_ind++)
	{
		unsigned prefixLength = 0;
	
		double graph_size = graph_dataset[g_ind].vertexCount + graph_dataset[g_ind].edgeCount; // size of graph g_ind

		{
			// For (g_ind)th graph prefix length will be this only for computation with any other graph.

			double invUbound = (double)1.0/ubound; // inverse of ubound
			prefixLength = 1 +(unsigned)(ceil((graph_size*(double)(1.0-invUbound)))); // Prefix Length
		}
		
		// Constructing the rank-list 
	
		// Putting vertex and edge ranks together
		vector<unsigned> graph_ranks;
		for(int vtx_ind = 0; vtx_ind < graph_dataset[g_ind].vertices.size(); vtx_ind++)
		{
			pair<unsigned, unsigned> vtx_pair = make_pair(graph_dataset[g_ind].vertices[vtx_ind], graph_dataset[g_ind].vertices[vtx_ind]);
			graph_ranks.push_back(rank[vtx_pair]);
		}
		
		for(int edge_ind = 0; edge_ind < graph_dataset[g_ind].edges.size(); edge_ind++)
		{
			pair<unsigned, unsigned> edge_pair;
			if(graph_dataset[g_ind].edges[edge_ind].first > graph_dataset[g_ind].edges[edge_ind].second)
				edge_pair = make_pair(graph_dataset[g_ind].edges[edge_ind].second, graph_dataset[g_ind].edges[edge_ind].first);
			else
				edge_pair = make_pair(graph_dataset[g_ind].edges[edge_ind].first, graph_dataset[g_ind].edges[edge_ind].second);
			graph_ranks.push_back(rank[edge_pair]);
		}

		// sort the ranks of the graph in descending order  
		sort(graph_ranks.begin(), graph_ranks.end(), greater <>());


		// crop graph's rank-list upto prefix-length
		for(int pref = 0; pref < graph_size; pref++)
			rankList[g_ind].push_back(graph_ranks[pref]);
		// rankList is till prefix length

		if(false)
		{
			bucket[g_ind].resize(no_of_buckets);
			for(int buck_ind = 0; buck_ind < no_of_buckets; buck_ind++)
				bucket[g_ind][buck_ind].resize(graph_size, 0);

			vector<unsigned> sumBucket(no_of_buckets, 0);

			// traversing graph's rank-list in ascending order 
			for(int grank_ind = graph_size-1; grank_ind >= 0; grank_ind--)
			{
				sumBucket[graph_ranks[grank_ind]%no_of_buckets]++;
				for(int buck_ind = 0; buck_ind < no_of_buckets; buck_ind++)
				{
					bucket[g_ind][buck_ind][grank_ind] = sumBucket[buck_ind];
				}
			}
		}
	}
*/
}
int temp=0;

void VEO:: calculate_sparse_table(vector<Graph> &graph_dataset, int g_ind, long double minPrevSize)  // INVERTED INDEXING FOR RANKLISTS
{
		chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

	double graph_size = graph_dataset[g_ind].vertexCount + graph_dataset[g_ind].edgeCount; 
	double invUbound = (double)1.0/ubound; // inverse of ubound
	int	prefixLength = 1 +(unsigned)(ceil((graph_size*(double)(1.0-invUbound)))); // Prefix Length
	sparse_table.clear(); 	// re-setting the sparse table

	///////////////////////////////////////////////////
	// for this particular graph g_ind, make a sparse table with help of inverted list

	for(int i = 0; i < prefixLength; i++){ // traversing a graph

		int rank = rankList[g_ind][i];
		
		for(int j = InvertedIndex[rank].size()-1; j>=0; j--){ // traversing a rank's inverted list(whisch is asc order) so we traverse from last

			int gr = InvertedIndex[rank][j].first;
			long double PrevSize = graph_dataset[gr].vertexCount + graph_dataset[gr].edgeCount;
			if(PrevSize < minPrevSize)		// loose filter condition to break and stop unnecessary computation
				break;

			if(gr == g_ind)		// skipping itself in inverted list
				continue;
			
			if( sparse_table.count( gr ) == 0)  // if occuring for the first time
				sparse_table[ gr ].first = 1;
			else
				{
					sparse_table[ gr ].first ++;		// incrementing the count of commons in 2 graphs list upto a length
					sparse_table[ gr ].second = InvertedIndex[rank][j].second;  // storing that length of gr (shorter graph) upto which commons have been counted
				}
			}
	}

	///////////////////////////////////////////////////

	//Now lets create InvertedIndex for (g_ind)th graph..
	
	for(int pref = 0; pref < prefixLength; pref++)
		InvertedIndex[ rankList[g_ind][pref] ].push_back({g_ind, pref+1}); // pushing graph_id and position of the particular attribute/rank in that graph's list 

chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
/*	temp +=(unsigned long long int)(1e-6*chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count());
	cout<<temp<<" ";*/
	//}
}

// index each input graphs in dataset
void VEO:: index(vector<Graph> &graph_dataset, int mode, bool isBucket, int no_of_buckets)
{	
	ranking(graph_dataset); // ranking vertices and edges together
	buildPrefix(graph_dataset, mode, isBucket, no_of_buckets);
}


// Applies prefix filter on Graph g1 and g2
bool VEO:: PrefixFilter(Graph &g1, Graph &g2, int index1, int index2, int mode, bool isBucket, int no_of_buckets, long unsigned &indexCount, double threshold)
{
	bool out = true;
	// atleast 1 common you got in the prefix part, then .
    if(sparse_table.count(index2) != 0)// or sparse_table[index2].count(index1) != 0)  // this sparse table is of g1/index1 graph
    {
		out=false;
	}

	if(out)
		return out;
	indexCount++;

	return out;
}


// Applies prefix filter on Graph g1 and g2
bool VEO:: PositioningFilter(Graph &g1, Graph &g2, int index1, int index2, int mode, bool isBucket, int no_of_buckets, long unsigned &indexCount, double threshold)
{
	unsigned size1 = g1.vertexCount + g1.edgeCount; // Size of Graph g1
	unsigned size2 = g2.vertexCount + g2.edgeCount; // Size of Graph g2

	unsigned prefix1;
	unsigned prefix2;

		double invUbound = (double)1.0/ubound; // inverse of ubound
		//prefix1 = rankList[index1].size(); // prefix-length of graph g1
		prefix1 = 1 +(unsigned)(ceil((size1*(double)(1.0-invUbound)))); // Prefix Length
		//prefix2 = rankList[index2].size(); // prefix-length of graph g2
		prefix2 = 1 +(unsigned)(ceil((size2*(double)(1.0-invUbound)))); // Prefix Length

	unsigned start1 = 0;
	unsigned start2 = 0;
	bool out = true;
	long double commonTotal=0;

	// atleast 1 common you got in the prefix part, then .
    if(sparse_table.count(index2) != 0)// or sparse_table[index2].count(index1) != 0)  // this sparse table is of g1/index1 graph
    {
		out=false;
        int partial_score, remaining1, remaining2;

        partial_score = sparse_table[index2].first;

        remaining1 = size1 - prefix1;  // big graph
        remaining2 = size2 - sparse_table[index2].second;   // sparse_table[index2].second is the position upto which partial_score is stored in sparse table
        
        long double sizeSum = (long double)(size1 + size2);
        long double Common = (long double)1.0*(partial_score + min(remaining1, remaining2));
        long double veoEstimate = (long double)(200.0*Common)/(sizeSum);

	//	if(g1.gid == 5 and g2.gid == 33)
       // cout<<index1<<" "<<index2<<" "<<" "<<veoEstimate<<"\n\n";
        if((long double)veoEstimate < (long double)threshold)
            out = true;
    }

	if(out)
		return out;
	indexCount++;

	return out;
}

void VEO:: Preprocess_Suffix(vector<Graph> &graph_dataset, int no_of_buckets)
{
	graph_bucket.resize(graph_dataset.size()+1);
	cout<<"Max rank is "<<r<<"\n";	
	cout<<"no_of_buckets is "<<no_of_buckets<<"\n";
	int bucket_size = ceil(1.0*r/no_of_buckets);
	cout<<"buckets size "<<bucket_size<<"\n";

	for(int g=0; g<graph_dataset.size(); g++){
		graph_bucket[g].resize(no_of_buckets+1,0);

		for(int r=0; r< rankList[g].size(); r++){
			graph_bucket[ g ][ ceil(1.0*rankList[g][r]/bucket_size) ]++;
		}
		/*if(g==8 or g==1)
		{cout<<"Buckets for "<<graph_dataset[g].gid<<" is :";
		for(auto j:graph_bucket[g])cout<<j<<" ";
		cout<<"\n";}*/
	}
}


bool VEO:: SuffixFilter(Graph &g1, Graph &g2, int index1, int index2, double threshold, bool isBucket, int no_of_buckets, long unsigned &SuffixFilterCount){

	unsigned size1 = g1.vertexCount + g1.edgeCount; // Size of Graph g1
	unsigned size2 = g2.vertexCount + g2.edgeCount; // Size of Graph g2
	double invUbound = (double)1.0/ubound; // inverse of ubound
	unsigned prefix1 = (unsigned)(ceil((size1*(double)(1.0-invUbound)))); // Prefix Length
	unsigned prefix2 = sparse_table[index2].second; // Prefix Length
	unsigned partial_score = sparse_table[index2].first;
	bool out = true;

	while(prefix1>=0 && prefix2>=1 && rankList[index1][prefix1] != rankList[index2][prefix2-1])
		prefix1--;
//	double H_max = size1 + size2 /*-((size1 + size2)*threshold/100.0) + 2*(sparse_table[index1][index2].first)*/ - 2.0*ceil((size1 + size2)*(threshold/200.0));// - (prefix1 + sparse_table[index1][index2].second);
//	cout<<"1\n";
//	double H = suffix_filter(rankList[index1], rankList[index2], prefix1, size1-1, sparse_table[index2].second, size2-1,H_max, 1);
//cout<<"2\n";
//cout<<H<<" "<<H_max<<" "<<size1 + size2-prefix1-sparse_table[index1][index2].second<<"\n";
/*	if(H <= H_max)
		return false;

	return true;*/

	int bucket_size = ceil(1.0*r/no_of_buckets);
	int current_bucket = ceil(1.0*rankList[index1][prefix1]/bucket_size); // current bucket of g1's prefix junction 
	int end_pt = (current_bucket - 1)*bucket_size + 1; 	// end pt of current bucket of g1's prefix junction
//cout<<index1<<" "<<index2<<" "<<size1<<" "<<prefix1<<" "<<rankList[index1][prefix1]<<" "<<current_bucket<<" "<<end_pt<<" "<<partial_score<<" ";
	
	/*if(g1.gid==27 and g2.gid==32){//cout<<"savdhan :"<<index1<<" "<<index2<<"\n";}
		cout<<"Ranklist for 27 is :";
		for(auto j:rankList[11])cout<<j<<" ";
		cout<<"\n";
		cout<<"Ranklist for 32 is :";
		for(auto j:rankList[0])cout<<j<<" ";
		cout<<"\n";
		cout<<"prefixes : "<<prefix1<<" "<<prefix2<<"\n";
		cout<<"end_pt : "<<end_pt<<"\n";
		cout<<"actual : "<<intersection_vertices(g1.vertices,g2.vertices)<<" "<<intersection_edges(g1.edges,g2.edges)<<" ";
		cout<<"partial : "<<partial_score<<" ";

		}*/
	
	
	
	// Linear traversal till nearest next bucket
	
	while(prefix1<rankList[index1].size() && prefix2<rankList[index2].size() && rankList[index1][prefix1] >= end_pt && rankList[index2][prefix2] >= end_pt)
	{
		if(rankList[index1][prefix1] == rankList[index2][prefix2])
		{
			partial_score++;
			prefix1++;
			prefix2++;
		}
		else
		if(rankList[index1][prefix1] > rankList[index2][prefix2])
			prefix1++;
		else
			prefix2++;
	}
//cout<<partial_score<<" ";
	// Now, calculating suffix overlap from buckets
	for(int i=current_bucket-1; i>=1; i--){
		partial_score += min(graph_bucket[index1][i], graph_bucket[index2][i]);
	}
//	cout<<partial_score<<" \n";
//unsigned actual =intersection_vertices(g1.vertices,g2.vertices)+intersection_edges(g1.edges,g2.edges);
	
	//if(partial_score < actual ){cout<<index1<<" "<<index2<<" "<<partial_score<<" "<<actual<<" "<<++ans<<"\n";}

		long double sizeSum = (long double)(g1.vertexCount + g1.edgeCount + g2.vertexCount + g2.edgeCount);
		long double veoEstimate =(long double)200.0*((long double)partial_score/(long double)(sizeSum));
		if(veoEstimate < (long double)threshold)
			out = true;
		else
			{
				out = false;
				SuffixFilterCount++;
			}
		
	return out;
}


// MisMatching Filter
bool VEO:: mismatchingFilter(Graph &g1, Graph &g2, double &common, double threshold)
{
	int index1 = 0;
	int index2 = 0;
	unsigned degreeSum12 = 0;
	unsigned degreeSum21 = 0;
	while(index1 < g1.vertexCount && index2 < g2.vertexCount)
	{
		if(g1.vertices[index1] == g2.vertices[index2])
		{
			common++;
			index1++;
			index2++;
		}
		else if(g1.vertices[index1] > g2.vertices[index2])
		{
			degreeSum21 += g2.degrees[g2.vid_to_ind[g2.vertices[index2]]];
			index2++;
		}
		else
		{
			degreeSum12 += g1.degrees[g1.vid_to_ind[g1.vertices[index1]]];
			index1++;
		}
	}

	double ES = (double)min((g1.edgeCount - degreeSum12), (g2.edgeCount - degreeSum21) );
	ES = max(0.0,ES);
	double simEstimate = (double)(200.0*(double)((common + ES)/(double)(g1.edgeCount + g1.vertexCount + g2.edgeCount + g2.vertexCount)));
	return (simEstimate <= threshold);
}


bool VEO:: VertexLabelFilter(Graph &g1, Graph &g2, int index1, int index2, unordered_set<unsigned>& v_label, unordered_set<unsigned>& e_label, long unsigned &indexcount, double threshold){

	double common = 0.0;
	unsigned size1 = g1.vertexCount + g1.edgeCount; // Size of Graph g1
	unsigned size2 = g2.vertexCount + g2.edgeCount; // Size of Graph g2
	bool out = false;

	for(auto label : v_label){
		common += min(g1.VertexLabelMap[label], g2.VertexLabelMap[label]);
	}

	common += min((int)g1.edges.size(),(int)g2.edges.size());
	
	long double sizeSum = (long double)(size1 + size2);
    long double veoEstimate = (long double)(200.0*common)/(sizeSum);

	if((long double)veoEstimate < (long double)threshold)
            out = true;

	if(!out)
		indexcount++;

	return out;
}

bool VEO:: EdgeLabelFilter(Graph &g1, Graph &g2, int index1, int index2, unordered_set<unsigned>& v_label, unordered_set<unsigned>& e_label, long unsigned &indexcount, double threshold){

	double common = 0.0;
	unsigned size1 = g1.vertexCount + g1.edgeCount; // Size of Graph g1
	unsigned size2 = g2.vertexCount + g2.edgeCount; // Size of Graph g2
	bool out = false;

	for(auto label : v_label){
		common += min(g1.VertexLabelMap[label], g2.VertexLabelMap[label]);
	}

	for(auto label : e_label){
		common += min(g1.EdgeLabelMap[label], g2.EdgeLabelMap[label]);
	}
	
	long double sizeSum = (long double)(size1 + size2);
    long double veoEstimate = (long double)(200.0*common)/(sizeSum);

	if((long double)veoEstimate < (long double)threshold)
            out = true;

	if(!out)
		indexcount++;
		
	return out;
}