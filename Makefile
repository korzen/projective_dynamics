CC       ?= gcc
CXX      ?= g++
CFLAGS   = -std=c11 -g -I./ -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native `pkg-config --cflags gtk+-3.0 epoxy`
CXXFLAGS = -std=c++14 -g -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native
CUDA_LIBS= -L/usr/local/cuda-7.5/lib64/ -lcudart
LDFLAGS  = -lm `pkg-config --libs gtk+-3.0 epoxy`
EIGEN3   = `pkg-config --cflags eigen3`
VIENNACL = -I/home/sci/will/Downloads/ViennaCL-1.7.1/bin/include/


all: build_dir libpd_solver.so
	#$(CXX) $(CXXFLAGS) $(EIGEN3) -c src/backend/pd_eigen.cpp -o obj/pd_solver.o
	$(CC) $(CFLAGS) -c src/pd_io.c -o obj/pd_io.o
	$(CC) $(CFLAGS) -c src/pd_linalg.c -o obj/pd_linalg.o
	$(CC) $(CFLAGS) -c src/main.c -o obj/main.o
	$(CXX) $(CFLAGS) obj/pd_io.o obj/pd_linalg.o obj/main.o -o pd $(LDFLAGS) -L. -lpd_solver

libpd_solver.so: src/backend/pd_viennacl.cpp
	nvcc -x cu -arch=sm_20 -std=c++11 -O3 -DVIENNACL_WITH_CUDA $(VIENNACL) \
		--shared -Xcompiler -fPIC src/backend/pd_viennacl.cpp -o libpd_solver.so $(CUDA_LIBS)

.PHONY: build_dir
build_dir:
	mkdir -p obj

.PHONY: clean
clean:
	rm -rf obj; rm -f pd; rm -f libpd_solver.so

