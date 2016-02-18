CC       = gcc
CFLAGS   = -std=c11 -g -I./ -O3 -Wall `pkg-config --cflags gtk+-3.0 epoxy`
CXXFLAGS = -std=c++14 -g -O3 -Wall


all:
	$(CXX) $(CXXFLAGS) -c src/pd_solver.cpp -o obj/pd_solver.o
	$(CC) $(CFLAGS) obj/pd_solver.o src/pd_io.c src/main.c -o pd `pkg-config --libs gtk+-3.0 epoxy`
