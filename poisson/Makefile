override CXXFLAGS += -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0

all: poisson.exe

poisson.exe: poisson_test.cc poisson.cc
	    ${MPICXX} ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} ${LDLIBS} -o poisson.exe poisson.cc poisson_test.cc

clean:
	rm -f *.o *.exe
