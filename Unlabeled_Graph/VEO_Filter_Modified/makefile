# Define required macros here
SHELL = /bin/sh

CFLAG = -Wall -g
CC = g++
INCLUDE =
LIBS = -lm

Size = 2195
val = 2

all:

	${CC} -std=c++14 filter.cpp veo.cpp graph.cpp -o filter

	$ for N in `seq 80 1 100` ; do \
		./filter ./dataset/data_${val}k.txt 1 $$N ${Size} stat_${val}k;  \
	done; \

	$ for N in `seq 80 1 100` ; do \
		./filter ./dataset/data_${val}k.txt 2 $$N 0 ${Size} stat_${val}k;  \
	done; \

	$ for N in `seq 80 1 100` ; do \
		./filter ./dataset/data_${val}k.txt 3 $$N 0 0 ${Size} stat_${val}k;  \
	done; \

	$ for N in `seq 80 1 100` ; do \
		./filter ./dataset/data_${val}k.txt 4 $$N 0 0 ${Size} stat_${val}k ; \
	done; \

	$ for N in `seq 80 1 100` ; do \
		./filter ./dataset/data_${val}k.txt 5 $$N 0 20 ${Size} stat_${val}k ; \
	done;

clean:
	rm -f *.o *.o
