# Dataset generation by adding NDG by modifying an existing data file. (without labels)

Here graphs are being modified by using modification.txt
For now, graph after being parsed are nowhere written back. We are just doing modifications in a graph randomly.
If Needed we can add these modified graphs (near duplicates) in actual graph data.

**Imp Obs:
Many times random function cannot generate a few numbers paricularly, reason (being random), thus causing program to go infinite. I have bounded such modifications to 4 times of number of vertices.

compile : g++ duplicatGraphsGeneration.cpp  -o GraphGeneration

run : ./GraphGeneration dataset.txt modifications.txt
