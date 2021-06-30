#include<bits/stdc++.h>
using namespace std;

int main(){


    ifstream inp("data_test.txt");
    ofstream otp("../data/data_test1.txt");
    
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

            int vertex_label;
            inp >> tag >> vertex_label;
            otp << tag << " " << i << " " << vertex_label <<"\n";

            mp[vertex_label] = i;
        }

        for(int i=0; i<edge_count; i++){

            int u, v;
            inp >> tag >> u >> v;
            u = mp[u];
            v = mp[v];

            otp << tag << " " << u << " " << v << " " << 25000 <<"\n";
        }
    }

    return 0;
}