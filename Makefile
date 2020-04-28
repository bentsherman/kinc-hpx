BINS = kinc-omp kinc-hpx

all: ${BINS}

%-omp: %-omp.cpp
	g++ -o $@ $^ -std=c++11 -fopenmp

%-hpx: %-hpx.cpp
	g++ -o $@ $^ \
		$(shell pkg-config --cflags --libs hpx_application_release) \
		-DHPX_APPLICATION_NAME=$@ \
		-lhpx_iostreams

clean:
	rm -f *.o ${BINS}
