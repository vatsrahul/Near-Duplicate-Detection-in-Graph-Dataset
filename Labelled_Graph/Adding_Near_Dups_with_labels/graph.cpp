#include "graph.h"

 // adds an edge to the graph
void Graph:: pushEdge(unsigned src_vtx,unsigned dest_vtx, int ec){		// vid are pushed
	if(src_vtx > dest_vtx)
		{
			this->edges.push_back(make_pair(src_vtx, dest_vtx));
			this->eid_to_ec[{src_vtx, dest_vtx}]=ec;
		}
	else {
		this->edges.push_back(make_pair(dest_vtx, src_vtx));
		this->eid_to_ec[{dest_vtx, src_vtx}]=ec;
		}
}

// reads the graph from input file 
void Graph:: readGraph(istream &inp){
	char tag;
	// first line should be of the format "g vertexCount(unsigned int) edgeCount(unsigned int) gid(unsigned int)"
	inp >> tag; // the tag 'g'
	inp >> vertexCount; // the no. of vertices in the graph
	inp >> edgeCount; // the no. of edges in the graph
	inp >> gid; // the graph-id of the graph
	vertices.resize(vertexCount);
	
	unsigned ind = 0;
	unsigned vid, src_vtx, dest_vtx, ec;
	char vc;  //Vertex type -- vc 
	for(int vtx_ind=0; vtx_ind < vertexCount; vtx_ind++)
	{
		// each line for each vertex should be in the format like: "v vid(unsigned int)"
		inp >> tag >> vid >> vc; // the tag 'v' along with the vertex-id
		vertices[vtx_ind] = vid;
		vid_to_ind[vid] = vtx_ind; // mapping vertex-id to its index
		vid_to_vc[vid] = vc;

		max_label = max(max_label, vc-'A');     // this is max label in a graph, used when we add a vertex
	}
	for(int e_ind = 0; e_ind < edgeCount; e_ind++)
	{
		// each line for each edge should be in the format like: "e vid_src(unsigned int) vid_dest(unsigned int)"
		inp >> tag >> src_vtx >> dest_vtx >> ec; // the tag 'e' along with the source and destination vertex-ids
		// Undirected graph : adding edge source to destination and destination to source
		pushEdge(src_vtx, dest_vtx, ec);			
	}
}

// prints details of the graph
/*void Graph::displayGraph()
{
	cout <<"g "<< gid << ":" << endl;
	cout <<"Vertex Count: "<< vertexCount << endl;
	cout <<"Edge Count: "<< edgeCount << endl;
	for(int vtx_ind = 0; vtx_ind<vertexCount; vtx_ind++)
	{
		cout << "v" << vtx_ind << ": " << vertices[vtx_ind] << endl;
	}
	for(int edg_ind = 0; edg_ind<edgeCount; edg_ind++)
	{
		cout << "e" << edg_ind << ": " << edges[edg_ind].first << " " << edges[edg_ind].second << endl;
	}
}*/