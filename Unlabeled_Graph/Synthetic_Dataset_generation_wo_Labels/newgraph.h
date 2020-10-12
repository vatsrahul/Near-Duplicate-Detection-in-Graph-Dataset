#pragma once
#include "common.h"

class Graph
{
	public:
		unsigned gid;
		vector<unsigned> vertices; // Vertex-set
		vector<pair<unsigned, unsigned> > edgeo; // Edge-Set 
		unordered_set<unsigned> useta;
		vector<pair<unsigned,unsigned>> edgeindex;
};

