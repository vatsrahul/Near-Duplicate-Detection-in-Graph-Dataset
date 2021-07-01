#include "veo.h"
#include <cmath>
#include "graph.h"
#include <bits/stdc++.h>



// Edge comparator
bool edgeComp(pair<unsigned, unsigned> &a, pair<unsigned, unsigned> &b,Graph &g1,Graph &g2)// since we have stored labels in edges, so we need not use vid_to_vc here.
{
	if(a.first == b.first)
		return a.second < b.second;
	return a.first < b.first;
}

// Intersection of vertex lists of 2 graphs
double intersection_vertices(vector<unsigned> &s1, vector<unsigned> &s2,unsigned gid1,unsigned gid2,Graph &g1,Graph &g2)
{


	unsigned s1_iter = 0;
	unsigned s2_iter = 0;
	double common = 0;

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
	}

	return common;
}

// Intersection of edge lists of 2 graphs
double intersection_edges(vector<unsigned> &s1, vector<unsigned> &s2,unsigned gid1,unsigned gid2,Graph &g1,Graph &g2)
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

// VEO Similarity computation
double computeSimilarity(Graph &g1, Graph &g2, double &commonV)
{
	if(commonV == 0)
		commonV = intersection_vertices(g1.vertices, g2.vertices,g1.gid,g2.gid,g1,g2);
	commonV+=intersection_edges(g1.edges, g2.edges,g1.gid,g2.gid,g1,g2);

	//cout << commonV << " "<<g1.vertexCount+g2.vertexCount+g1.edgeCount+g2.edgeCount<<"\n";

	double simScore = (double)(200.0*(commonV)/(double)(g1.vertexCount+g2.vertexCount+g1.edgeCount+g2.edgeCount));

	return simScore;
}
