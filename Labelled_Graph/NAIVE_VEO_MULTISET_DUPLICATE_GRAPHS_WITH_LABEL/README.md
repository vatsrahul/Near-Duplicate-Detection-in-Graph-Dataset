# Naive VEO on Duplicate graph Dataset with Labels

**Changes** i am making according to my current uderstanding.
1.  Vertices array have ids in them which are sorted, and then in finding intersection labels are used.
	example: 1,C 1,H 2,O 2,C 3,N 3,O........

	But it should be sorted acc to labels and not vertex ids.
	like : 1,C 2,C 1,H 3,N 2,O 3,O......
	Then only we can apply merge 2 sort list algo to find intersection correctly.

2.  For Edges though, labels are stored in pair format (large label, small label) in edges array.
	and then intersection is applied but vid_to_vc is used. (A BLUNDER), i removed vid_to_vc as already labels were pushed in edges.



Program Flow:

	1. Parse the input file for Graph Dataset using parseGraphDataset().

	2. Sorts the graph dataset:
	
	3. For each graph pair among all nC2 pairs:
		i) compute the similarity: computeSimilarity() function in veo.cpp
			simScore = 100*{2*[commonV/(g1.vertexCount+g2.vertexCount+g1.edgeCount+g2.edgeCount)]}

	4. Record statistics like similarity-score frequency, similarity-time for each graph.

	5. Store all graph pairs which have similarity score greater than given threshold in a file.

	6. Record and store final statistics like total similarity time taken and total memory used for the whole dataset.

	7. Record frequency distribution of similarity score of graph pairs in freq_distr_file.txt.


Compiling this code:
$ g++ naive.cpp veo.cpp graph.cpp -o naive

Running the code:
$ ./naive graph-dataset-file simScore_threshold dataset-size res-dir

For Example:
$ ./naive ./dataset/data_file.txt 60 10 stat

Input Arguements:

	1. graph-dataset-file:
		Requirements:
			1. dataset must be already sorted based on LABELS.
			2. file should be in the following format:
			
				g vertexCount edgeCount graph_id1
				v vertex_id1 label
				v vertex_id2 label
				.
				.
				.
				e source_vertex_id1 label destination_vertex_id1 label
				e source_vertex_id2 label destination_vertex_id2 label
				.
				.
				.
				g vertexCount edgeCount graph_id2
				.
				.
				.

	2. simScore_threshold: only those graph pairs with similarity score greater than this threshold will be stored in a file

	3. dataset-size: size of dataset to be selected

	4. result-directory: directory in which all the output files will be stored

Output Files:

	The program creates a directory named 'res-dir' in which it has 2 files:

		1. all_graph_file.txt which has all graph pairs with similarity greater than simScore_threshold as follows:
			
			graph_id1 graph_id2 similarity_score

		2. stat_final.txt which will have all that stats as given below:

			GSimJoin: VEO Similarity(naive)
			Dataset size: 1000
			Similarity Score Threshold: 60
			Similar Graph Pairs: 190
			Memory used: 3.81641 MB
			Similarity Time: 0 milliseconds
			Total Time Taken: 264 milliseconds

		3. graph_details: it stores all the statistics related to a particular graph(each file of the form: "g_graph-index_graph-id_sim.txt"):
			g_graph-index_graph-id_sim.txt contents:
				no. of graphs with simScore greater than the threshold
				graph-id graph-id2 simScore
				.
				.
				.
				followed by 102 lines of simScore frequency for that graph
				similarity time for that graph
	
		4. freq_distr_file.txt

			102 lines of pair of nos as follows
			<sim-score> <sim-pair-count>
			where,
				the first line indicates no of graph pairs with similarity score 0
				the last line indicates no of isomorphic graph pairs 
				and the rest of the line indicates no of graph pairs with similarity score sim-score to sim-score+1

