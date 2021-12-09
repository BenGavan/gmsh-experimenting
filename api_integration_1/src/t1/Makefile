#
#output: main.o message.o
#	g++ main.o message.o -o output

main.o: main.cpp
	g++ -std=c++11 -o main main.cpp -lgmsh

#g++ -std=c++11 -o main main.cpp -lgmsh

run: main.cpp
	g++ -std=c++11 -o main main.cpp -lgmsh
	./main

clean:
	rm *.o output