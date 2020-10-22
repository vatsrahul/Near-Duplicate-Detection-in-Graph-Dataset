#include "common.h"
using namespace std;

class Vertex{

	public:
		unsigned vid; // Vertex id
		double quality; // Vertex Quality (PageRank)
		vector<Vertex*> edges; // List of Edges adjacent to the vertex along with its quality values
		unsigned label;

		Vertex()
		{
			vid = 0;
			quality = 0;
			label = 0;
		}
};


class Graph{

	public:
		unsigned gid; // graph-id 
		unsigned vertexCount; // No. of Vertices
		unsigned edgeCount; // No. of Edges
		vector<Vertex> vertices;
		unordered_map<unsigned,unsigned> vid_to_ind; // Mapping between Vertex Label and its Index in adjacancy list 
		unordered_map<unsigned,char> vid_to_label;
		Graph(){
			gid = 0;
			vertexCount = 0;
			edgeCount = 0;
		}
		void displayGraph(); // displays the graph's details
		void readGraph(ifstream &inp); // Parsing an input file to create the graph
		void sortGraph(); // Sort a graph vertices: based on quality of vertices and edges: based on quality of destn vts

};

