CC       ?= gcc
CXX      ?= g++
CFLAGS   = -std=c11 -g -I./ -O3 -Wall -march=native `pkg-config --cflags gtk+-3.0 epoxy`
CXXFLAGS = -std=c++14 -g -O3 -Wall -march=native
LDLIBS   = -lm


all:
	$(CC) $(CFLAGS) -c src/pd_solver.c -o obj/pd_solver.o
	$(CC) $(CFLAGS) obj/pd_solver.o src/pd_io.c src/pd_linalg.c src/main.c -o pd `pkg-config --libs gtk+-3.0 epoxy` $(LDLIBS)
