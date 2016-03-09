CC       ?= gcc
CXX      ?= g++
CFLAGS   = -std=c11 -g -I./ -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native `pkg-config --cflags gtk+-3.0 epoxy`
CXXFLAGS = -std=c++14 -g -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native
LDFLAGS  = -lm `pkg-config --libs gtk+-3.0 epoxy`
EIGEN3   = `pkg-config --cflags eigen3`


all: build_dir pd
	$(CXX) $(CXXFLAGS) $(EIGEN3) -c src/backend/pd_eigen.cpp -o obj/pd_solver.o
	$(CC) $(CFLAGS) -c src/pd_io.c -o obj/pd_io.o
	$(CC) $(CFLAGS) -c src/pd_linalg.c -o obj/pd_linalg.o
	$(CC) $(CFLAGS) -c src/main.c -o obj/main.o
	$(CXX) $(CFLAGS) obj/pd_solver.o obj/pd_io.o obj/pd_linalg.o obj/main.o -o pd $(LDFLAGS)

.PHONY: build_dir
build_dir:
	mkdir -p obj

.PHONY: clean
clean:
	rm -rf obj; rm pd

