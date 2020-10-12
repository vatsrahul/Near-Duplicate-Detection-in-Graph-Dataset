# Dataset generation including labels without near duplicates

This is same as approach.cpp of Synthetic Dataset generation + here vertices have an associated label with it just like chemical compounds. 


compile: g++ generation.cpp -o control
run : ./control 10 500 400 80 42 26  data_1.txt

(In general, ./control graph_count max_vertices min_verices max_%_edges min_%_edges label_count file_name