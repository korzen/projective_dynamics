CC       ?= gcc
CXX      ?= g++
CFLAGS   = -std=c11 -g -I./ -O3 -Wall -Wextra -fstrict-aliasing -march=native
CXXFLAGS = -std=c++14 -g -I./ -Iext/ -O3 -Wall -Wextra -fopenmp -fstrict-aliasing -march=native `pkg-config --cflags epoxy glfw3`
CUDA_LIBS= -L/usr/local/cuda-7.5/lib64/ -lcudart -lcusolver -lcusparse
LDFLAGS  = -lm `pkg-config --libs glfw3 epoxy` -l:ext/jsmn/libjsmn.a
EIGEN3   = `pkg-config --cflags eigen3`

ifeq ($(BACKEND), cuda)
	LIBPD_SOLVER=pd_solver_cuda
	PD_SO_NAME=libpd_solver_cuda.so
	SOLVER_DEFINES=-DVIENNACL_WITH_CUDA
else ifeq ($(BACKEND), cusparse_high)
	LIBPD_SOLVER=pd_solver_cuda
	PD_SO_NAME=libpd_solver_cuda.so
	SOLVER_DEFINES=-DVIENNACL_WITH_CUDA -DUSE_CUSPARSE=1
else ifeq ($(BACKEND), cusparse_low)
	LIBPD_SOLVER=pd_solver_cuda
	PD_SO_NAME=libpd_solver_cuda.so
	SOLVER_DEFINES=-DVIENNACL_WITH_CUDA -DUSE_CUSPARSE=1 -DUSE_CUSPARSE_LOW_LEVEL=1
else ifeq ($(BACKEND), opencl)
	LIBPD_SOLVER=pd_solver_opencl
	PD_SO_NAME=libpd_solver_opencl.so
else
	LIBPD_SOLVER=pd_solver_eigen
	PD_SO_NAME=libpd_solver_eigen.so
endif

all: build_dir $(PD_SO_NAME) build_jsmn
	$(CXX) $(CXXFLAGS) -c src/pd_io.c -o obj/pd_io.o
	$(CXX) $(CXXFLAGS) -c src/pd_linalg.c -o obj/pd_linalg.o
	$(CXX) $(CXXFLAGS) -c src/main.c -o obj/main.o
	$(CXX) $(CXXFLAGS) -Iext/imgui -c ext/imgui/imgui_impl_glfw_gl3.cpp -o obj/imgui_impl_glfw_gl3.o
	$(CXX) $(CXXFLAGS) -c ext/imgui/imgui.cpp -o obj/imgui.o
	$(CXX) $(CXXFLAGS) -c ext/imgui/imgui_draw.cpp -o obj/imgui_draw.o
	$(CXX) $(CXXFLAGS) obj/pd_io.o obj/pd_linalg.o obj/main.o obj/imgui_impl_glfw_gl3.o obj/imgui.o obj/imgui_draw.o -o pd $(LDFLAGS) -L. -l $(LIBPD_SOLVER)

libpd_solver_eigen.so: src/backend/pd_eigen.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN3) --shared -fPIC $^ -o $@

libpd_solver_cuda.so: src/backend/pd_viennacl.cpp
	nvcc -x cu -arch=compute_52 -code=sm_52 -std=c++11 -O3 $(SOLVER_DEFINES) \
		--shared -Xcompiler -fPIC src/backend/pd_viennacl.cpp -o $@

libpd_solver_opencl.so: src/backend/pd_viennacl.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN3) --shared -fPIC -DVIENNACL_WITH_OPENCL $^ -o $@

build_jsmn:
	cd ext/jsmn && $(MAKE) -f Makefile

.PHONY: build_dir
build_dir:
	mkdir -p obj

.PHONY: clean
clean:
	rm -rf obj; rm -f pd; rm -f *.so

