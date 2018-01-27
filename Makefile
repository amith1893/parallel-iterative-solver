CXX=mpicxx
CXX1=c++11

SOURCE=jacobi.cc poisson.cc

all:jacobi.exe

%.exe:%.cc
	${CXX} -I. -std=${CXX1} -g ${SOURCE} -o $@

clean:
	rm -rf *.exe
