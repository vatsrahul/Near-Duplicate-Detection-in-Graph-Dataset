#include<bits/stdc++.h>
using namespace std;

int main(){


    ifstream inp("data_2k.txt");
    ofstream otp("../data/data_2k.txt");
    
    int size;
    inp >> size;

    while(size--){

        int vertex_count, edge_count, graph_name;
        char tag;
        unordered_map<int,int> mp;
        mp.clear();

        inp >> tag >> vertex_count >> edge_count >> graph_name;

        otp<<"t # "<<graph_name<<"\n";

        for(int i=0; i<vertex_count; i++){

            int vertex_id;
            char vertex_label;
            inp >> tag >> vertex_id >> vertex_label;
            otp << tag << " " << i << " " << vertex_label-'A' <<"\n";

            mp[vertex_id] = i;
        }

        for(int i=0; i<edge_count; i++){

            int u, v, edge_label;
            inp >> tag >> u >> v >> edge_label;
            u = mp[u];
            v = mp[v];

            otp << tag << " " << u << " " << v << " " << edge_label <<"\n";
        }
    }

    return 0;
}