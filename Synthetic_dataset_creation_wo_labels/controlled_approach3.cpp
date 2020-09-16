#include<iostream>
#include<vector>
#include<unordered_map>
#include<unordered_set>
#include<time.h>
#include<fstream>
#include<bits/stdc++.h>
#include "graph.h"


using namespace std;

//One way of generating near duplicates is that keep track of uset of edges considered for geneating graph like this you have for each graph a uset now what we do is that vary like last two elements of the present array try to generate the graps for these and introduce them into dataset but this requires lot of space  
//Then we will store edge index,uset as a data member of class graph 


void createGraph(unsigned gnumber,unsigned &nv,unsigned &ne,unsigned long maxV,unsigned long  minV,unsigned long maxE,unsigned long minE,vector<unsigned>& vertex_ids,unordered_map<unsigned,unordered_set<unsigned> >& edges,vector<pair<unsigned,unsigned>> edgeindex,unordered_set<unsigned> useta, int r);

//Generating Random Values 
unsigned getRandomVertexCount(unsigned maxV,unsigned minV);
unsigned getRandomVertexLabel(unsigned maxV);
unsigned getRandomPercentageofedges(unsigned maxE,unsigned minE);
unsigned getRandomEdgeLabel(unsigned maxE);

//Printing and Writing to resultant files
void printGraph(unsigned gnumber,unsigned nv,unsigned  ne,vector<unsigned>& vertex_ids,unordered_map<unsigned,unordered_set<unsigned> >& edges);
void writeToFile(unsigned gnumber,unsigned nv,ofstream &datafile,vector<unsigned>& vertex_ids,unordered_map<unsigned,unordered_set<unsigned> >& edges);
void writeToGraphs(int i,vector<unsigned> &vertex_ids,unordered_map<unsigned,unordered_set<unsigned> > &edges,vector<Graph> &graph_dataset);

//Generating  vertix ids for all vertices in present graph  
void getVertex_ids(unsigned  nv,vector<unsigned> &vertex_ids);




int main(int argc, char const *argv[])
{
    
   //Input Constrains 
    int NGraphs = stoi(argv[1]);
    if(NGraphs <= 0)
    {
        cerr<<"Input required: valid number of graphs"<<endl;  //To validate desired Inputs 
		exit(0);
    }
    
    unsigned long maxV =   stoul(argv[2]);
    unsigned long  minV =   stoul(argv[3]);
    unsigned long maxE =   stoul(argv[4]);
    unsigned long minE =   stoul(argv[5]);
    
	 
    srand((unsigned)time(NULL));

    unsigned nv,ne; //No of Vertices and Edges in the current Graph 
    ofstream datafile(argv[6]);
    datafile<<NGraphs<<"\n"; 		// writing no f graphs in file

    vector<Graph> graph_dataset(NGraphs);  //Initializing Graph Dataset
 
    for(int g = 0; g < NGraphs; g++)
   {
	     vector<unsigned> vertex_ids; //VectorOfVectorLabels
	     unordered_map<unsigned,unordered_set<unsigned> > edges;//unordered_map of final edges of graph
	     vector<pair<unsigned,unsigned>> edgeindex;//To store indices of the correspoing edges 
	     unordered_set<unsigned> useta;   // utility to ensure same edge is not selected out of nv*(nv-1) edges
		 createGraph(g+1,nv,ne,maxV,minV,maxE,minE,vertex_ids,edges,edgeindex,useta,g);//Creaters Graph with given Constrains
	     writeToFile(g+1,nv,datafile,vertex_ids,edges);//Writing graph to resultant dataset file
	     writeToGraphs(g,vertex_ids,edges,graph_dataset);//writing to the Graph Datastrucure
   }
	datafile.close();
}


unsigned getRandomVertexCount(unsigned long maxV,unsigned long minV){
	return rand()%(maxV-minV)+minV;
}
unsigned getRandomVertexLabel(unsigned long maxV){
	return rand()%(maxV)+1;
}
unsigned getRandomEdgeLabel(unsigned long edgecount){
	return rand()%(edgecount)+0;
}
unsigned getRandomPercentageofedges(unsigned long maxE,unsigned long minE){
	return rand()%(maxE-minE)+minE;
}
bool isEdgeExist(unsigned i,unsigned j,unordered_map<unsigned,unordered_set<unsigned> > &edges){
	return edges.find(i)!=edges.end() && edges[i].find(j)!=edges[i].end();
}
/*Returns true if edge is already exist otherwise returns false*/


void getVertex_ids(unsigned nv,vector<unsigned> &vertex_ids,unsigned long maxV){

	unordered_set<unsigned> uset; //To Check that no vid is repeated 
	unsigned vid;
	while(nv){
		vid  = getRandomVertexLabel(maxV);  //we should get random vertices they need not be inside nv nv is "Total No" of vertices
		if(uset.find(vid)==uset.end()){
			uset.insert(vid);
			vertex_ids.push_back(vid);
			nv--;			
		}
	}
	
}

void createGraph(unsigned gnumber,unsigned &nv,unsigned &ne,unsigned long maxV,unsigned long minV,unsigned long maxE,unsigned long minE,vector<unsigned>& vertex_ids,unordered_map<unsigned,unordered_set<unsigned> >& edges,vector<pair<unsigned,unsigned>> edgeindex,unordered_set<unsigned> useta,int r){
    
  
    nv = getRandomVertexCount(maxV,minV); //To get a random vertex count between min no vertices and max no of vertices
    getVertex_ids(nv,vertex_ids,maxV);
    unsigned NEdges = nv*(nv-1); //No of total no of edges possible for current number of vertices 

    unsigned currEP = getRandomPercentageofedges(maxE,minE); //The percentage of edges that are allowded to form the graph  
    
    unsigned curr_no_ed = ((currEP * NEdges)/100);//The percentage of edges that are allowded to form the graph

    
	unsigned long ecount = 0;
	for(unsigned i=0;i<nv && ecount < NEdges;i++)
	{
		for(unsigned j=0;j<nv && ecount< NEdges; j++)
		{
                     if(i != j)
		     {
				edgeindex.push_back({vertex_ids[i],vertex_ids[j]}); // storing all possible edges except self loops
				ecount++;
            }
		}
	}
	//cout<<"edges "<<ecount<<"\n";
	
	while(curr_no_ed > 0)
	{ 
		//cout << curr_no_ed << '\n';
		int e = getRandomEdgeLabel(ecount);
	    if(useta.find(e) == useta.end()){
		    
		    useta.insert(e);

            if(!isEdgeExist(edgeindex[e].second, edgeindex[e].first, edges))
			{
		    	    edges[edgeindex[e].first].insert(edgeindex[e].second);

				 curr_no_ed--;
			}
		}
	}
}

void writeToFile(unsigned gnumber,unsigned nv,ofstream &datafile,vector<unsigned>& vertex_ids,unordered_map<unsigned,unordered_set<unsigned> >& edges){
	
	unsigned ne = 0;
	for(auto mitr=edges.begin();mitr!=edges.end();mitr++)
	for(auto sitr=mitr->second.begin();sitr!=mitr->second.end();sitr++)
		ne++;

	datafile << "g "<<nv<<" "<<ne<<" "<<gnumber<<endl;

	for(auto vid:vertex_ids)
		datafile << "v " << vid << endl;
	for(auto mitr=edges.begin();mitr!=edges.end();mitr++){
		for(auto sitr=mitr->second.begin();sitr!=mitr->second.end();sitr++){
			datafile << "e " << mitr->first <<" "<< *sitr <<endl;
		}
	}
}

void printGraph(unsigned gnumber,unsigned nv,unsigned ne,vector<unsigned>& vertex_ids,unordered_map<unsigned,unordered_set<unsigned> >& edges){ //printing the graph 
	cout << "g "<<gnumber << " vertices: " << nv << " edges: " << ne <<endl;
	for(auto vid:vertex_ids)
		cout << "v "<<vid<<endl;
	for(auto mitr=edges.begin();mitr!=edges.end();mitr++){
		for(auto sitr=mitr->second.begin();sitr!=mitr->second.end();sitr++){
			cout<<mitr->first<<" "<<*sitr<<endl;
		}
	}
}


void writeToGraphs(int i,vector<unsigned> &vertex_ids,unordered_map<unsigned,unordered_set<unsigned> > &edges,vector<Graph> &graph_dataset)//Writing it to Graph
{
	

	for(auto vid:vertex_ids)
	{
        	graph_dataset[i].vertices.push_back(vid);
        }
	
        for(auto mitr=edges.begin();mitr!=edges.end();mitr++){
		for(auto sitr=mitr->second.begin();sitr!=mitr->second.end();sitr++){
			graph_dataset[i].edgeo.push_back({mitr->first,*sitr});
		}
	}
	/*for(auto v:graph_dataset[i].vertices)
		cout<< v << " ";
	cout<<endl;

	for(auto e:graph_dataset[i].edgeo){
		
		cout << e.first << " " << e.second;
                cout << '\n';		
	}
	cout<<endl;
        cout << endl; */
  
}


