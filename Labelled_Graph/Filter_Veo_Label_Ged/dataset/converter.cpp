#include<bits/stdc++.h>
using namespace std;
 
int main(){


    ifstream inp("aids_dataset.txt");  // add a 'q' at end of the file temporarily.
    ofstream otp("data_aids.txt");    // after running add size at top

        unsigned size=0;
        char tag;

        inp >> tag;

        while(tag == 't' and tag != '\n'){
        
            unsigned long i=0,vertex_count, edge_count, graph_name;
            vector<unsigned> vertex_id(2000), vertex_label(2000), u(2000), v(2000), edge_label(2000);
            size++;

            inp >> tag >> graph_name; 
            inp >> tag;
            while(tag == 'v'){
                inp >> vertex_id[i] >> vertex_label[i];
                inp >> tag;
                i++;
            }
            vertex_count = i;

            i=0;
            while(tag == 'e'){
                inp >> u[i] >> v[i] >> edge_label[i];
                inp >> tag;
                i++;
            }
            edge_count = i;

            otp<< "g " << vertex_count <<" "<< edge_count << " " <<graph_name<<"\n";
            cout<< "g " << vertex_count <<" "<< edge_count << " " <<graph_name<<"\n";
            
            
            for(int i=0; i<vertex_count; i++){
                otp << "v " << vertex_id[i] <<" "<< vertex_label[i] <<"\n";
            }
             
            for(int i=0; i<edge_count; i++){
                otp << "e " << u[i] <<" "<< v[i] <<" "<< edge_label[i] <<"\n";
            }
        
        }

        cout<<"Graph size "<<size <<"\n Add it above dataset\n";

    return 0;
}