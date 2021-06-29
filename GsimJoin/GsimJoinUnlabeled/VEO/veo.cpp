// Run using the following command
// g++ veo.cpp

#include <bits/stdc++.h>

using namespace std;

#define DATASET 5500

int intersect(vector<int> a,vector<int> b){
    int n = a.size();
    int m = b.size();

    int i=0,j=0;
    int count = 0;
    while(i<n && j<m){
        if(a[i]<b[j]) i++;
        else if(a[i]>b[j]) j++;
        else{
            count++;
            i++;
            j++;
        }
    }
    return count;
}

int main(){
    ifstream fp;
    int t;
    char c;
    int u,v,e,idx;
    vector<vector<int>> gVertex(DATASET+1);
    vector<vector<int>> gEdges(DATASET+1); 
    fp.open("../data/data_2k.txt");
    while(fp>>c){
        switch(c){
            case 'g': fp>>idx>>u>>t;
                      if(t>1){
                          sort(gVertex[t-1].begin(),gVertex[t-1].end());
                          sort(gEdges[t-1].begin(),gEdges[t-1].end());
                      }
                      break;
            case 'v': fp>>v;
                      gVertex[t].push_back(v);
                      break;
            case 'e': fp>>u>>v;
                      if(u< v){
                          idx = (u-1)*100007 + v-1;
                      }
                      else{
                          idx = (v-1)*100007 + u-1;
                      }
                      gEdges[t].push_back(idx);
                      break;
        }
    }
    fp.close();
    int n = gVertex.size();
    sort(gVertex[DATASET].begin(),gVertex[DATASET].end());
    sort(gVertex[DATASET].begin(),gVertex[DATASET].end());
    

    ofstream ofp("../Output/result.txt");
    ifstream ifp("../run/2k_result.txt");
    ofp<<fixed<<setprecision(3);
int i,j,ged,count=0;
 long double mx=0.0, mn = 1.0,sum=0.0;
vector<long double> vec_ged[6];
          
    while(ifp>>i and  ifp>>j and ifp>>ged){
		count++;
            int tot = gVertex[i].size() + gVertex[j].size() + gEdges[i].size() + gEdges[j].size();
            int vertexIntsct = intersect(gVertex[i],gVertex[j]);
            int edgeIntsct = intersect(gEdges[i],gEdges[j]);
            long double sim = (long double)(2)*(long double)(vertexIntsct+edgeIntsct); 
	    if(tot > 0){
              sim = sim/tot;
            }
            else{
              sim = 0.0;
            }
            mx = max(mx,sim);
            mn = min(mn,sim);
            sum += sim;
                ofp<<i<<" "<<j<<" "<<ged<<" ";
                ofp<<sim<<endl;
		vec_ged[ged].push_back(sim);
            
        }
ofp<<"\nMax sim : "<<mx<<endl;
ofp<<"Min sim : "<<mn<<endl;
ofp<<"avg sim : "<<1.0*sum/count<<endl<<endl;

for(int i=1;i<=5;i++){
	for(auto j: vec_ged[i])
		ofp<<j<<" ";
ofp<<endl;
}

  ifp.close();
    ofp.close();
    cout<<endl;   
    return 0;
}
