#include "poisson.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <sys/time.h>

#define pi 3.141592653589793

int main(int argc, char* argv[]) {
  std::vector<point_t> x;
  MPI_Init(&argc, &argv);

  int r;
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  if (argc < 2) {
    if (r == 0) printf("Usage: %s n\n", argv[0]);
    MPI_Finalize();
    exit(1);
  }
  int n = atoi(argv[1]);

  MPI_Comm grid_comm;
  MPI_Comm start_comm = MPI_COMM_WORLD;
  poisson_setup(start_comm, n, grid_comm, x);
  std::vector<real_t> a(0);
  std::vector<real_t> u(0);
  std::vector<real_t> lu(0);
  for (auto p : x) {
    a.push_back(12);
    u.push_back(sin(4*pi*p.x) * sin(10*pi*p.y) * sin(14*pi*p.z));
  }
  lu.resize(u.size());
 
  struct timeval tic, toc;
  gettimeofday(&tic, NULL);
  laplacian(grid_comm, n, a, u, lu);
  gettimeofday(&toc, NULL);
  double elapsed_time = (toc.tv_sec - tic.tv_sec) + (toc.tv_usec - tic.tv_usec)*1e-6;
  MPI_Finalize();

  if (r==0) printf("Time (s): %e\n", elapsed_time);
}

