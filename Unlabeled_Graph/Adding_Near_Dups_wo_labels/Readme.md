# Dataset generation by adding NDG by modifying an existing data file. (without labels)

Here graphs are being modified by using modification.txt
For now, graph after being parsed are nowhere written back. We are just doing modifications in a graph randomly.
If Needed we can add these modified graphs (near duplicates) in actual graph data.

**Imp Obs:
Many times random function cannot generate a few numbers paricularly, reason (being random), thus causing program to run for long. I have bounded such modifications to 4 times of number of max vertices.**
After discussio with sir, came to know this will happen when >90% edges,  but we will takr around 50-60% edges. So no problem.

Above problem solved by linearly choosing vertex. For edges i still prefer not to do any edit. (as it wud be n^2).

**Max vertex taken = 30**

**From data file some graphs randomly will be modified (based on modifications file) and then written in dupicates file. Graph id for duplicates are generated sequentially further. So u need to just copy duplicates.txt in last of data.txt after running below 2 commands. REMEBER TO CAHNGE TOTAL NO OF GRAPHS AT TOP.**

modifications.txt for :

2k
10
20 1
20 2
20 3
20 4
20 5
20 6
20 7
20 8
20 9
20 10

5k
10
40 1
40 2
40 3
40 4
40 5
30 6
30 7
30 8
30 9
30 10

10k
10
50 1
50 2
50 3
50 4
50 5
50 6
50 7
50 8
50 9
50 10



**compile** : g++ duplicatGraphsGeneration.cpp  -o GraphGeneration

**run** : ./GraphGeneration dataset.txt modifications.txt duplicates.txt
