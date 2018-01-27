#include "poisson.h"
#include <iostream>
#include <cmath>
#include <cassert>
#include <sys/time.h>

void get_rank_and_size(const MPI_Comm &comm, int &n, int &r) {
  MPI_Comm_size(comm, &n);
  MPI_Comm_rank(comm, &r);
}

void get_my_coord_bounds(const MPI_Comm &grid_comm, int r, int n, int coord_starts[],
    int coord_ends[]) {
  int coords[DIM] = {0};
  int dims[DIM] = {0};
  int trash[DIM] = {0};
  MPI_Cart_get(grid_comm, DIM, dims, trash, coords);
  MPI_Cart_coords(grid_comm, r, DIM, coords);
  for (int i=0; i<DIM; i++) {
    coord_starts[i] = n*coords[i]/dims[i];
    coord_ends[i] = n*(1+coords[i])/dims[i];
  }
}

void poisson_setup(MPI_Comm &comm, int n, MPI_Comm &grid_comm,
    std::vector<point_t> &x) {
  int dims[DIM] = {0};
  int periods[DIM] = {0};
  int n_p;
  int r;

  get_rank_and_size(comm, n_p, r);

  MPI_Dims_create(n_p, DIM, dims);
  MPI_Cart_create(comm, DIM, dims, periods, 0, &grid_comm);

  int coord_starts[DIM];
  int coord_ends[DIM];
  get_my_coord_bounds(grid_comm, r, n, coord_starts, coord_ends);

  point_t next_pt;
  x.clear();
  for (int i=coord_starts[0]; i<coord_ends[0]; i++) {
    for (int j=coord_starts[1]; j<coord_ends[1]; j++) {
      for (int k=coord_starts[2]; k<coord_ends[2]; k++) {
        next_pt.x = static_cast<real_t>(i)/(n-1);
        next_pt.y = static_cast<real_t>(j)/(n-1);
        next_pt.z = static_cast<real_t>(k)/(n-1);
        x.push_back(next_pt);
      }
    }
  }
}

//void poisson_residual(MPI_Comm &grid_comm, int n, std::vector<real_t> &a,
    //std::vector<real_t> &f, std::vector<real_t> &v, real_t &res)
void laplacian(const MPI_Comm &grid_comm, int n, const std::vector<real_t> &a,
                   const std::vector<real_t> &u, std::vector<real_t> &lu) {
  std::vector<real_t> r_h(u.size(), 0);
  std::vector<real_t> send_buffers[DIM][2];
  std::vector<real_t> recv_buffers[DIM][2];
  int neighbors[DIM][2];
  int coord_starts[DIM];
  int coord_ends[DIM];

  int r;
  int n_p;
  get_rank_and_size(grid_comm, n_p, r);
  get_my_coord_bounds(grid_comm, r, n, coord_starts, coord_ends);

  int mi = coord_ends[0] - coord_starts[0];
  int mj = coord_ends[1] - coord_starts[1];
  int mk = coord_ends[2] - coord_starts[2];

  assert(mi*mj*mk == u.size());

  real_t h = 1.0/n;

#if 0
  for (int i=0; i<mi*mj*mk; i++) {
	  //lu[i] = (a[i]*u[i] - r_h[i]);
	  printf("u[i] %lf\n", u[i]);
  }
  printf("==============================================================================================\n");
#endif
  // neighbor[0][0] has coords < coord_starts[0]  (i)
  // neighbor[0][1] has coords > coord_ends[0]
  // neighbor[1][0] has coords < coord_starts[1]  (j)
  // neighbor[1][1] has coords > coord_ends[1]
  // neighbor[2][0] has coords < coord_starts[2]  (k)
  // neighbor[2][1] has coords > coord_ends[2]
  for (int i=0; i<DIM; i++) {
    MPI_Cart_shift(grid_comm, i, 1, &(neighbors[i][0]), &(neighbors[i][1]));
  }

  for (int i=0; i<mi; i++) {
    for (int j=0; j<mj; j++) {
      send_buffers[2][0].push_back(u[    i*mj*mk +      j*mk +      0]);
      send_buffers[2][1].push_back(u[    i*mj*mk +      j*mk + (mk-1)]);
    }
  }
  for (int i=0; i<mi; i++) {
    for (int k=0; k<mk; k++) {
      send_buffers[1][0].push_back(u[    i*mj*mk +      0*mk +      k]);
      send_buffers[1][1].push_back(u[    i*mj*mk + (mj-1)*mk +      k]);
    }
  }
  for (int j=0; j<mj; j++) {
    for (int k=0; k<mk; k++) {
      send_buffers[0][0].push_back(u[     0*mj*mk +     j*mk +      k]);
      send_buffers[0][1].push_back(u[(mi-1)*mj*mk +     j*mk +      k]);
    }
  }

  MPI_Request is_req[6];
  MPI_Request ir_req[6];
  for (int dim=0; dim<DIM; dim++) {
    for (int dir=0; dir<2; dir++) {
      size_t s = send_buffers[dim][dir].size();
      recv_buffers[dim][dir].assign(s, 0);
      MPI_Isend(send_buffers[dim][dir].data(), s, MPI_DOUBLE,
          neighbors[dim][dir], 0, grid_comm, &(is_req[2*dim+dir]));
      MPI_Irecv(recv_buffers[dim][dir].data(), s, MPI_DOUBLE,
          neighbors[dim][dir], 0, grid_comm, &(ir_req[2*dim+dir]));
    }
  }

  for (int i=0; i<mi; i++) {
    for (int j=0; j<mj; j++) {
      for (int k=0; k<mk; k++) {
        r_h[i*mj*mk + j*mk + k] = -6*u[i*mj*mk + j*mk + k];
        if (i>0) {
          r_h[i*mj*mk + j*mk + k] += u[(i-1)*mj*mk + j*mk + k];
        }
        if (i<mi-1) {
          r_h[i*mj*mk + j*mk + k] += u[(i+1)*mj*mk + j*mk + k];
        }
        if (j>0) {
          r_h[i*mj*mk + j*mk + k] += u[i*mj*mk + (j-1)*mk + k];
        }
        if (j<mj-1) {
          r_h[i*mj*mk + j*mk + k] += u[i*mj*mk + (j+1)*mk + k];
        }
        if (k>0) {
          r_h[i*mj*mk + j*mk + k] += u[i*mj*mk + j*mk + (k-1)];
        }
        if (k<mj-1) {
          r_h[i*mj*mk + j*mk + k] += u[i*mj*mk + j*mk + (k+1)];
        }
      }
    }
  }

  MPI_Status ir_stat[6];
  MPI_Waitall(6, ir_req, ir_stat);

  for (int i=0; i<mi; i++) {
    for (int j=0; j<mj; j++) {
      r_h[i*mj*mk + j*mk + 0] += recv_buffers[2][0][i*mj+j];
      r_h[i*mj*mk + j*mk + (mk-1)] += recv_buffers[2][1][i*mj+j];
    }
  }
  for (int i=0; i<mi; i++) {
    for (int k=0; k<mk; k++) {
      r_h[i*mj*mk + 0*mk + k] += recv_buffers[1][0][i*mk+k];
      r_h[i*mj*mk + (mj-1)*mk + k] += recv_buffers[1][1][i*mk+k];
    }
  }
  for (int j=0; j<mj; j++) {
    for (int k=0; k<mk; k++) {
      r_h[0*mj*mk + j*mk + k] += recv_buffers[0][0][j*mk+k];
      r_h[(mi-1)*mj*mk + j*mk + k] += recv_buffers[0][1][j*mk+k];
    }
  }

  for (int i=0; i<mi*mj*mk; i++) {
    lu[i] = (a[i]*u[i] - r_h[i]);
	//printf("lu[i] %lf\n", lu[i]);
  }
}
