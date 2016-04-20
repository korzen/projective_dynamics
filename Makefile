CC       ?= gcc
CXX      ?= g++
CFLAGS   = -std=c11 -O3 -I./ -Wall -Wextra -fstrict-aliasing -march=native
CXXFLAGS = -std=c++14 -O3 -I./ -Iext/ -Wall -Wextra -fopenmp -fstrict-aliasing \
		   -march=native `pkg-config --cflags epoxy glfw3`
CUDA     ?= /usr/local/cuda/
LDFLAGS  = -lm `pkg-config --libs glfw3 epoxy` -l:ext/jsmn/libjsmn.a
EIGEN3   = `pkg-config --cflags eigen3`
BACKENDF = active_backend.txt

ifeq ($(BACKEND), cuda)
	LIBPD_SOLVER=pd_solver_cuda
	PD_SO_NAME=libpd_solver_cuda.so
	SOLVER_DEFINES=-DVIENNACL_WITH_CUDA
	CUDA_LIBS=-L$(CUDA)/lib64/ -lcudart
else ifeq ($(BACKEND), cusparse_high)
	LIBPD_SOLVER=pd_solver_cuda
	PD_SO_NAME=libpd_solver_cuda.so
	SOLVER_DEFINES=-DVIENNACL_WITH_CUDA -DUSE_CUSPARSE=1
	CUDA_LIBS=-L$(CUDA)/lib64/ -lcudart -lcusolver -lcusparse
else ifeq ($(BACKEND), cusparse_low)
	LIBPD_SOLVER=pd_solver_cuda
	PD_SO_NAME=libpd_solver_cuda.so
	SOLVER_DEFINES=-DVIENNACL_WITH_CUDA -DUSE_CUSPARSE=1 -DUSE_CUSPARSE_LOW_LEVEL=1
	CUDA_LIBS=-L$(CUDA)/lib64/ -lcudart -lcusolver -lcusparse
else ifeq ($(BACKEND), opencl)
	LIBPD_SOLVER=pd_solver_opencl
	PD_SO_NAME=libpd_solver_opencl.so
else
	BACKEND=eigen
	LIBPD_SOLVER=pd_solver_eigen
	PD_SO_NAME=libpd_solver_eigen.so
endif

# See if the backend has changed from what we built previously and force a rebuild/link
ifneq ($(BACKEND), $(shell test -f $(BACKENDF) && cat $(BACKENDF) || echo "$(BACKEND)" > $(BACKENDF)))
	UPDATE_BACKEND=$(file > $(BACKENDF),$(BACKEND))
endif

all: $(UPDATE_BACKEND) build_dir build_jsmn pd pd_benchmark pd_mesh_binary

pd: obj/pd_io.o obj/pd_linalg.o obj/main.o obj/libimgui.so $(PD_SO_NAME) $(BACKENDF)
	$(CXX) $(CXXFLAGS) obj/pd_io.o obj/pd_linalg.o obj/main.o -o pd $(LDFLAGS) \
	    -Wl,-rpath,./obj -L./obj/ -l imgui -L. -l $(LIBPD_SOLVER) $(CUDA_LIBS)

pd_benchmark: obj/benchmark.o obj/pd_io.o obj/pd_linalg.o $(PD_SO_NAME) $(BACKENDF)
	$(CXX) $(CXXFLAGS) obj/pd_io.o obj/pd_linalg.o obj/benchmark.o -o $@ $(LDFLAGS) \
		-Wl,-rpath,./obj -L./obj/ -l imgui -L. -l $(LIBPD_SOLVER) $(CUDA_LIBS)

pd_mesh_binary: obj/mesh_binary.o obj/pd_io.o
	$(CXX) $(CXXFLAGS) $^ -o $@ -L./ext/jsmn -ljsmn

obj/libimgui.so: obj/imgui/imgui_impl_glfw_gl3.o obj/imgui/imgui.o obj/imgui/imgui_draw.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^

obj/%.o: src/%.c
	$(CXX) $(CXXFLAGS) -c $^ -o $@

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

obj/imgui/%.o: ext/imgui/%.cpp
	$(CXX) -Iext/imgui -fPIC -c $^ -o $@

libpd_solver_eigen.so: src/backend/pd_eigen.cpp obj/libimgui.so $(BACKENDF)
	$(CXX) $(CXXFLAGS) $(EIGEN3) --shared -fPIC $< -o $@

libpd_solver_cuda.so: src/backend/pd_viennacl.cpp obj/libimgui.so $(BACKENDF)
	nvcc -x cu -std=c++11 -O3 $(SOLVER_DEFINES) \
		--shared -Xcompiler -fPIC -Iext/ $< -o $@ -L./obj/ -limgui

libpd_solver_opencl.so: src/backend/pd_viennacl.cpp $(BACKENDF)
	$(CXX) $(CXXFLAGS) $(EIGEN3) --shared -fPIC -DVIENNACL_WITH_OPENCL $< -o $@

build_jsmn:
	cd ext/jsmn && $(MAKE) -f Makefile

# Create the build dirs and an empty active backend file if
# it doesn't exist
.PHONY: build_dir
build_dir:
	mkdir -p obj; mkdir -p obj/imgui

.PHONY: clean
clean:
	rm -rf obj; rm -f pd; rm -f *.so; rm -f $(BACKENDF)

