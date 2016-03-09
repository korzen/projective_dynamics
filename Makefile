CC       ?= gcc
CXX      ?= g++
CFLAGS   = -std=c11 -g -I./ -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native `pkg-config --cflags gtk+-3.0 epoxy`
CXXFLAGS = -std=c++14 -g -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native
CUDA_LIBS= -L/usr/local/cuda-7.5/lib64/ -lcudart
LDFLAGS  = -lm `pkg-config --libs gtk+-3.0 epoxy`
EIGEN3   = `pkg-config --cflags eigen3`
# TODO Set VIENNACL through your environment or such
VIENNACL ?= -I/home/sci/will/Downloads/ViennaCL-1.7.1/bin/include/

ifeq ($(BACKEND), cuda)
	LIBPD_SOLVER=pd_solver_cuda
else ifeq ($(BACKEND), opencl)
	LIBPD_SOLVER=pd_solver_opencl
else 
	LIBPD_SOLVER=pd_solver_eigen
endif

all: build_dir $(LIBPD_SOLVER)
	$(CC) $(CFLAGS) -c src/pd_io.c -o obj/pd_io.o
	$(CC) $(CFLAGS) -c src/pd_linalg.c -o obj/pd_linalg.o
	$(CC) $(CFLAGS) -c src/main.c -o obj/main.o
	$(CXX) $(CFLAGS) obj/pd_io.o obj/pd_linalg.o obj/main.o -o pd $(LDFLAGS) -L. -l $(LIBPD_SOLVER)

pd_solver_eigen: src/backend/pd_eigen.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN3) --shared -fPIC $^ -o libpd_solver_eigen.so

pd_solver_cuda: src/backend/pd_viennacl.cpp
	nvcc -x cu -arch=sm_20 -std=c++11 -O3 -DVIENNACL_WITH_CUDA $(VIENNACL) \
		--shared -Xcompiler -fPIC src/backend/pd_viennacl.cpp -o libpd_solver_cuda.so $(CUDA_LIBS)

pd_solver_opencl: src/backend/pd_viennacl.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN3) --shared -fPIC -DVIENNACL_WITH_OPENCL $(VIENNACL) $^ -o libpd_solver_opencl.so

.PHONY: build_dir
build_dir:
	mkdir -p obj

.PHONY: clean
clean:
	rm -rf obj; rm -f pd; rm -f *.so

