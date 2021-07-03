#pragma once
#include "common.h"


class Graph
{
	public:
		unsigned gid, vertexCount, edgeCount; // graph-id, no. of vertices, no. of edges
		vector<unsigned> vertices; // Vertex-set
		vector<unsigned> degrees; // Degree-set
		unordered_map<unsigned,unsigned> vid_to_ind; // vid to index in adjacency list of graph
		unordered_map<unsigned,unsigned> vid_to_vc;
		vector< unsigned> edges; // Edge-Set 
		unordered_map<unsigned, unsigned> VertexLabelMap; // count of a type of vertex label
		unordered_map<unsigned, unsigned> EdgeLabelMap;  // count of a type of vertex label


		Graph(){
			gid = 0;
			vertexCount = 0;
			edgeCount = 0;
		}

		void readGraph(istream &inp, unordered_set<unsigned>& v_label, unordered_set<unsigned>& e_label); // reads the graph from input file 
		void pushEdge(unsigned u, unsigned v, unsigned ec, unordered_set<unsigned>& e_label); // adds an edge to the graph 
		void displayGraph(); // prints details of the graph

};

