# Define required macros here
SHELL = /bin/sh

CFLAG = -Wall -g
CC = g++
INCLUDE =
LIBS = -lm
all:
	${CC} -std=c++14 filter.cpp veo.cpp graph.cpp -o filter
	./filter ./dataset/data_5k.txt 1 60 5324 stat_5k
	./filter ./dataset/data_5k.txt 1 70 5324 stat_5k
	./filter ./dataset/data_5k.txt 1 80 5324 stat_5k
	./filter ./dataset/data_5k.txt 1 90 5324 stat_5k
	./filter ./dataset/data_5k.txt 1 95 5324 stat_5k
	./filter ./dataset/data_5k.txt 1 100 5324 stat_5k

	./filter ./dataset/data_5k.txt 2 60 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 2 70 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 2 80 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 2 90 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 2 95 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 2 100 0 5324 stat_5k

	./filter ./dataset/data_5k.txt 3 60 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 3 70 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 3 80 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 3 90 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 3 95 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 3 100 0 0 5324 stat_5k

	./filter ./dataset/data_5k.txt 4 60 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 4 70 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 4 80 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 4 90 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 4 95 0 0 5324 stat_5k
	./filter ./dataset/data_5k.txt 4 100 0 0 5324 stat_5k

	./filter ./dataset/data_5k.txt 5 60 0 20 5324 stat_5k
	./filter ./dataset/data_5k.txt 5 70 0 20 5324 stat_5k
	./filter ./dataset/data_5k.txt 5 80 0 20 5324 stat_5k
	./filter ./dataset/data_5k.txt 5 90 0 20 5324 stat_5k
	./filter ./dataset/data_5k.txt 5 95 0 20 5324 stat_5k
	./filter ./dataset/data_5k.txt 5 100 0 20 5324 stat_5k

clean:
	rm -f *.o *.o
