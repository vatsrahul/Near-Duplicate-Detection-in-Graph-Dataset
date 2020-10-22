#include "common.h"
#include "graph.h"
using namespace std;

bool vertexComp1(Vertex &v1, Vertex &v2)
{
	if(v1.quality == v2.quality)
		return v1.vid < v2.vid;
	return v1.quality > v2.quality;
}
bool vertexComp2(Vertex &v1, Vertex &v2)   // using only this comparator
{
	return v1.label < v2.label;    // sort on the basis of label
}

bool edgeComp(Vertex* &e1, Vertex* &e2) 
{
	if(e1->quality == e2->quality)
		return e1->vid < e2->vid;
	return e1->quality < e2->quality;
}


void Graph:: readGraph(ifstream &inp)
{
	char tag;
	// first line should be of the format "g vertexCount(unsigned int) edgeCount(unsigned int) gid(unsigned int)"
	inp >> tag; // first tag should be 'g'
	inp >> vertexCount;
	inp >> edgeCount;
	inp >> gid;

	vertices.resize(vertexCount);
	unsigned vid;
	char vc,sc,dc;
	for(int vtx_ind = 0; vtx_ind<vertexCount; vtx_ind++)
	{
		// each line for each vertex should be in the format like: "v vid(unsigned int)"
		inp >> tag >> vid >> vc;
		vertices[vtx_ind].vid = vid;	// vertex id
		vertices[vtx_ind].label = vc-'A'; // Label
		vid_to_label[vid]=vc;
	}
	sort(vertices.begin(),vertices.end(),vertexComp2);
	for(unsigned ind=0; ind<vertexCount; ind++)
	{
		vid_to_ind[vertices[ind].vid] = ind;			// Dont use this now.. anywhere...
	}
	int src_vtx, dest_vtx;

	for(int edg_ind = 0; edg_ind<edgeCount; edg_ind++)
	{
		// each line for each edge should be in the format like: "e vid_src(unsigned int) vid_dest(unsigned int)"
		inp >> tag >> src_vtx >> sc >> dest_vtx >> dc;
		edges.push_back(make_pair(src_vtx,dest_vtx));
		//vertices[vid_to_ind[src_vtx]].edges.push_back(&vertices[vid_to_ind[dest_vtx]]);
		//vertices[vid_to_ind[dest_vtx]].edges.push_back(&vertices[vid_to_ind[src_vtx]]);
	}
}


void Graph::displayGraph()	// NOT CHANGED FOR LABELS
{
	cout <<"g "<< gid << ":" << endl;
	cout <<"Vertex Count: "<< vertexCount << endl;
	cout <<"Edge Count: "<< edgeCount << endl;
	cout <<"graph shingle size: "<< shingles.size() << endl;
	cout <<"graph walk size: "<< walk.size() << endl;

	for(int vtx_ind = 0; vtx_ind<vertexCount; vtx_ind++)
	{
		cout << "v" << vtx_ind << ": " << vertices[vtx_ind].vid << endl;
		int edgeCount = vertices[vtx_ind].edges.size();
		for(int edg_ind = 0; edg_ind<edgeCount; edg_ind++)
		{
			cout << "e" << edg_ind << ": " << vertices[vtx_ind].vid << ", " << vertices[vtx_ind].edges[edg_ind]->vid << endl;
		}
	}
}

void Graph::sortGraph()		// NOT USED SO NOT CHANGED FOR LABELS
{
	// First sorted the vertices based on the vertex quality
	sort(vertices.begin(),vertices.end(),vertexComp1);
	for(unsigned ind=0; ind<vertexCount; ind++)
	{
		vid_to_ind[vertices[ind].vid] = ind;
	}
	// Then sort the edge list of each vertex based on destn vertex's quality
	for(unsigned i=0; i<vertexCount; i++)
	{
		sort(vertices[i].edges.begin(),vertices[i].edges.end(),edgeComp);
	}

}
