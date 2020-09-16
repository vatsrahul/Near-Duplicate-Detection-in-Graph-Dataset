Following repository aims to generate some random graph datasets (without labels).


controlled.cpp - It is a naive way of implementation where we are selecting edges based on 2 vertex_ids generated each time randomly. It seems efficient.
 
controlled_approach3 - It is a modified way where we generate all the edges (quadratic order) and select an edge at random using edge index and including it in the graph. (It is inefficient for a graph of say size a million)

**Interesting observation : 
1. In approach 3 whenever u want say 95% edges then, random function might take huge time to find an edge which is not taken before. so while using this approach keep in mind that. (50% edges is fine)

2. To generate large graph use controlled.cpp

3. In future multiple cores can be used to reduce runtime.



To Generate the Results Run the Following Commands - 

1.  g++ controlled.cpp(or corresponding file) -o control

2.  ./control NumberofGraphstoBeGenerated maximumnoofvertexices minimumnoofverices maxedgepercentage minedgepercentage datafile
    
     For Example ./control 10 1000 800 60 42 data_1.txt
