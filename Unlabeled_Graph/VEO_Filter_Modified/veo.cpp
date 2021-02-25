#include "veo.h"
#include <cmath>
#include <chrono>

bool freqComp(pair<pair<unsigned, unsigned>,unsigned> v1, pair<pair<unsigned, unsigned>,unsigned> v2)
{ 
    return (v1.second > v2.second); 
}

// Edge comparator
bool edgeComp(pair<unsigned, unsigned> &a, pair<unsigned, unsigned> &b)
{
	if(a.first == b.first)
		return a.second < b.second;
	return a.first < b.first;
}

// Intersection of vertex lists of 2 graphs
double intersection_vertices(vector<unsigned> &s1, vector<unsigned> &s2)
{
	unsigned s1_iter = 0;
	unsigned s2_iter = 0;
	double common = 0;

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
	return common;
}

// Intersection of edge lists of 2 graphs
double intersection_edges(vector<pair<unsigned, unsigned>> &s1, vector<pair<unsigned, unsigned>> &s2)
{
	unsigned s1_iter = 0;
	unsigned s2_iter = 0;
	double common = 0;

	for(auto s1_iter = s1.begin(), s2_iter = s2.begin(); s1_iter != s1.end() && s2_iter != s2.end(); )
	{
		if(edgeComp(*s1_iter,*s2_iter))
			s1_iter++;
		else if(edgeComp(*s2_iter,*s1_iter))
			s2_iter++;
		else
		{
			
			s1_iter++;
			s2_iter++;
			common++;		
		}
	}
	return common;
}

// VEO Similarity computation
double VEO:: computeSimilarity(Graph &g1, Graph &g2, double &commonV)
{
	if(commonV == 0)
		commonV = intersection_vertices(g1.vertices, g2.vertices);
	commonV+=intersection_edges(g1.edges, g2.edges);
	
	double simScore = (double)(200.0*(commonV)/(double)(g1.vertexCount+g2.vertexCount+g1.edgeCount+g2.edgeCount));
	return simScore;
}

// Sorts vertex and edge set of graph dataset
void VEO:: sortVertexEdge(vector<Graph> &graph_dataset)
{
	for(int i = 0; i < graph_dataset.size(); i++)
	{
		sort(graph_dataset[i].vertices.begin(), graph_dataset[i].vertices.end());
		sort(graph_dataset[i].edges.begin(), graph_dataset[i].edges.end());
	}
}

void VEO:: printlist(int k)
{
	for(int i = 0; i < rankList[k].size(); i++)
		cout << i  << ": "<< rankList[k][i] << endl;
}

//globaly ranks vertices and edges together
void VEO:: ranking(vector<Graph> &graph_dataset) 
{
	// traversing the graph-dataset
	for(int g_ind = 0; g_ind < graph_dataset.size(); g_ind++)
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
	unsigned long r = 1;
	for(auto entry = freqList.begin(); entry!=freqList.end(); entry++)
	{
		rank[entry->first] = r;
		r++;
	}

	// setting the size of invertecindex list as max value of rank - r
	InvertedIndex.resize(r);

}

void VEO:: buildPrefix(vector<Graph> &graph_dataset, int mode, bool isBucket, int no_of_buckets)
{
	rankList.resize(graph_dataset.size());
	if(isBucket)
		bucket.resize(graph_dataset.size());
	
	for(int g_ind = 0; g_ind < graph_dataset.size(); g_ind++)
	{
		unsigned prefixLength = 0;
	
		double graph_size = graph_dataset[g_ind].vertexCount + graph_dataset[g_ind].edgeCount; // size of graph g_ind

		if(mode == 3)	 // static Mode
		{
			// For (g_ind)th graph prefix length will be this only for computation with any other graph.

			double invUbound = (double)1.0/ubound; // inverse of ubound
			prefixLength = 1 +(unsigned)(ceil((graph_size*(double)(1.0-invUbound)))); // Prefix Length
		}
		else		// Dynamic Mode
		{
			// since it is dynamic, For (g_ind)th graph, prefix length will be calculated depending on other
			//  graphs size. so for now we will compute rank list for full size.

			prefixLength = graph_size;
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
		for(int pref = 0; pref < prefixLength; pref++)
			rankList[g_ind].push_back(graph_ranks[pref]);

		// rankList is till prefix length

		if(isBucket)
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


// INVERTED INDEXING FOR RANKLISTS

	sparse_table.resize(graph_dataset.size());

	for(int g_ind = 0; g_ind < graph_dataset.size(); g_ind++)
	{

		int prefixLength = rankList[g_ind].size();

		///////////////////////////////////////////////////
		// for this particular graph, make a sparse table with help of inverted list

	//	set<int> ss; // ss set contains list of graphs which are similar to (g_ind)th graph

		for(int i = 0; i < prefixLength; i++){ // traversing a graph

			int rank = rankList[g_ind][i];

			for(int j = 0; j < InvertedIndex[rank].size(); j++){ // traversing a rank's inverted list

				int gr = InvertedIndex[rank][j].first;
				if(gr == g_ind)		// skipping itself in inverted list
					continue;
				
				//ss.insert ( InvertedIndex[rank][j] );
				if( sparse_table[g_ind].count( gr ) == 0)  // if occuring for the first time
					sparse_table[g_ind][ gr ].first = 1;
				else
					{
						sparse_table[g_ind][ gr ].first ++;		// incrementing the count of commons in 2 graphs list upto a length
						sparse_table[g_ind][ gr ].second = InvertedIndex[rank][j].second;  // storing that length of gr (shorter graph) upto which commons have been counted
					}
				}
		}

		///////////////////////////////////////////////////

		//Now lets create InvertedIndex for (g_ind)th graph..
		
		for(int pref = 0; pref < prefixLength; pref++)
			InvertedIndex[ rankList[g_ind][pref] ].push_back({g_ind, pref+1}); // pushing graph_id and position of the particular attribute/rank in that graph's list 

	}

	// for(auto x: rankList[38])
	// cout<<x<<" ";
	// cout<<"\n";
	// for(auto x: rankList[26])
	// cout<<x<<" ";
	// cout<<"\n";
	

	/*cout<<"5 33 "<< sparse_table[38][26]<<"\n";
	cout<<intersection_vertices(graph_dataset[38].vertices,graph_dataset[26].vertices ) <<" - "<< 
	intersection_edges(graph_dataset[38].edges,graph_dataset[26].edges) <<"\n";
	cout<<"rlist :"<<rankList[38].size()<<"\n";
	cout<<"rlist :"<<rankList[26].size()<<"\n";
	
	cout<<"remain :"<<graph_dataset[38].edgeCount+graph_dataset[38].vertexCount - rankList[38].size()<<"\n";
	cout<<"remain :"<<graph_dataset[26].edgeCount+graph_dataset[26].vertexCount - rankList[26].size()<<"\n";
	

	cout<<"6 8 "<< sparse_table[8][6]<<"\n";
	cout<<intersection_vertices(graph_dataset[6].vertices,graph_dataset[8].vertices ) + 
	intersection_edges(graph_dataset[6].edges,graph_dataset[8].edges) <<"\n";*/

}

// index each input graphs in dataset
void VEO:: index(vector<Graph> &graph_dataset, int mode, bool isBucket, int no_of_buckets)
{	
	ranking(graph_dataset); // ranking vertices and edges together
	buildPrefix(graph_dataset, mode, isBucket, no_of_buckets);
}

// Applies index filter on Graph g1 and g2
bool VEO:: indexFilter(Graph &g1, Graph &g2, int index1, int index2, int mode, bool isBucket, int no_of_buckets, long unsigned &indexCount, long unsigned &partitionCount, double threshold)
{
	unsigned size1 = g1.vertexCount + g1.edgeCount; // Size of Graph g1
	unsigned size2 = g2.vertexCount + g2.edgeCount; // Size of Graph g2

	unsigned prefix1;
	unsigned prefix2;

	if(mode == 3) // static Mode
	{
		prefix1 = rankList[index1].size(); // prefix-length of graph g1
		prefix2 = rankList[index2].size(); // prefix-length of graph g2
	}
	if(mode == 4)	// Dynamic Mode
	{
		unsigned common = (unsigned)floor(((double)(threshold/200.0))*(double)(size1+size2));
		prefix1 = size1 - common + 1; // prefix-length of graph g1
		prefix2 = size2 - common + 1; // prefix-length of graph g2
	}

	unsigned start1 = 0;
	unsigned start2 = 0;
	bool out = true;
	long double commonTotal=0;

	// atleast 1 common you got in the prefix part, then .
    if(sparse_table[index1].count(index2) != 0 or sparse_table[index2].count(index1) != 0)
        //{out = false; cout<<"f\n";}
    {
		out=false;
        int partial_score, remaining1, remaining2;

        if(sparse_table[index1].count(index2) != 0)
            partial_score = sparse_table[index1][index2].first;
        else
            partial_score = sparse_table[index2][index1].first;

//cout<<partial_score<<"\n";
        remaining1 = size1 - prefix1;  // big graph
        remaining2 = size2 - sparse_table[index1][index2].second;   // sparse_table[index1][index2].second is the position upto which partial_score is stored in sparse table
        
        long double sizeSum = (long double)(size1 + size2);
        long double Common = (long double)1.0*(partial_score + min(remaining1, remaining2));
        long double veoEstimate = (long double)(200.0*Common)/(sizeSum);

	//	if(g1.gid == 5 and g2.gid == 33)
       // cout<<index1<<" "<<index2<<" "<<" "<<veoEstimate<<"\n\n";
        if((long double)veoEstimate < (long double)threshold)
            out = true;
    }


	/*while(start1 < prefix1 && start2 < prefix2)
	{
		if(rankList[index1][start1] == rankList[index2][start2] && isBucket)
		{
			out = false;
			start1++;
			start2++;
			commonTotal++;		
		}
		else if(rankList[index1][start1] == rankList[index2][start2] && !isBucket)
		{
			out = false;
			break;
		}
		else if(rankList[index1][start1] > rankList[index2][start2])
			start1++;
		else
			start2++;
	}*/
	if(out)
		return out;
	indexCount++;
	//cout << "prefix1: " << prefix1 << " , prefix2: " << prefix2 << endl;
	//cout << "start1: " << start1 << " , start2: " << start2 << endl;
	if(isBucket)
	{
		for(int i = 0; i < no_of_buckets; i++)
		{
			//cout << bucket[index1][i][start1] << " , " << bucket[index2][i][start2] << endl;
			commonTotal += (long double)min(bucket[index1][i][start1], bucket[index2][i][start2]);
		}

		long double sizeSum = (long double)(g1.vertexCount + g1.edgeCount + g2.vertexCount + g2.edgeCount);
		long double veoEstimate =(long double)200.0*((long double)commonTotal/(long double)(sizeSum));
		if(veoEstimate < (long double)threshold)
			out = true;
		else
		{
			out = false;
			partitionCount++;
		}
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
