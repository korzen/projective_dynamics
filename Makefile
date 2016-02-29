CC       ?= gcc
CXX      ?= g++
CFLAGS   = -std=c11 -g -I./ -O3 -Wall -Wextra -fstrict-aliasing -march=native `pkg-config --cflags gtk+-3.0 epoxy`
CXXFLAGS = -std=c++14 -g -O3 -Wall -Wextra -fstrict-aliasing -march=native
LDLIBS   = -lm


all:
	mkdir -p obj
#	$(CC) $(CFLAGS) -c src/backed/pd_verlet.c -o obj/pd_solver.o
	$(CXX) $(CXXFLAGS) `pkg-config --cflags eigen3` -c src/backend/pd_eigen.cpp -o obj/pd_solver.o
	$(CC) $(CFLAGS) -c src/pd_io.c -o obj/pd_io.o
	$(CC) $(CFLAGS) -c src/pd_linalg.c -o obj/pd_linalg.o
	$(CC) $(CFLAGS) -c src/main.c -o obj/main.o
	$(CXX) $(CFLAGS) obj/pd_solver.o obj/pd_io.o obj/pd_linalg.o obj/main.o -o pd `pkg-config --libs gtk+-3.0 epoxy` $(LDLIBS)
