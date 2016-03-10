CC       ?= gcc
CXX      ?= g++
CFLAGS   = -std=c11 -g -I./ -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native `pkg-config --cflags gtk+-3.0 epoxy`
CXXFLAGS = -std=c++14 -g -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native
CUDA_LIBS= -L/usr/local/cuda-7.5/lib64/ -lcudart
LDFLAGS  = -lm `pkg-config --libs gtk+-3.0 epoxy`
EIGEN3   = `pkg-config --cflags eigen3`

ifeq ($(BACKEND), cuda)
	LIBPD_SOLVER=pd_solver_cuda
	PD_SO_NAME=libpd_solver_cuda.so
else ifeq ($(BACKEND), opencl)
	LIBPD_SOLVER=pd_solver_opencl
	PD_SO_NAME=libpd_solver_opencl.so
else 
	LIBPD_SOLVER=pd_solver_eigen
	PD_SO_NAME=libpd_solver_eigen.so
endif

all: build_dir $(PD_SO_NAME)
	$(CC) $(CFLAGS) -c src/pd_io.c -o obj/pd_io.o
	$(CC) $(CFLAGS) -c src/pd_linalg.c -o obj/pd_linalg.o
	$(CC) $(CFLAGS) -c src/main.c -o obj/main.o
	$(CXX) $(CFLAGS) obj/pd_io.o obj/pd_linalg.o obj/main.o -o pd $(LDFLAGS) -L. -l $(LIBPD_SOLVER)

libpd_solver_eigen.so: src/backend/pd_eigen.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN3) --shared -fPIC $^ -o $@

libpd_solver_cuda.so: src/backend/pd_viennacl.cpp
	nvcc -x cu -arch=sm_20 -std=c++11 -O3 -DVIENNACL_WITH_CUDA \
		--shared -Xcompiler -fPIC src/backend/pd_viennacl.cpp -o $@ $(CUDA_LIBS)

libpd_solver_opencl.so: src/backend/pd_viennacl.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN3) --shared -fPIC -DVIENNACL_WITH_OPENCL $^ -o $@

.PHONY: build_dir
build_dir:
	mkdir -p obj

.PHONY: clean
clean:
	rm -rf obj; rm -f pd; rm -f *.so

