/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 */

#include <mpi.h>
#include <vector>

#define DIM 3
typedef double real_t;

union point_t
{
    real_t coord[3];
    struct {
        real_t x;
        real_t y;
        real_t z;
    };
};

/*
 * Set up the grid communicator for the Poisson problem and samples
 * the correct portion of the domain
 *
 * @comm input communicator, typically MPI_COMM_WORLD
 * @n problem size, there are n sample points along each axis
 * @grid_comm the grid communicator (output)
 * @x sample points assigned to the current MPI process based on its grid coordinates (output)
 *
 */
void poisson_setup(MPI_Comm &comm, int n, MPI_Comm &grid_comm, std::vector<point_t> &x);

/*
 * @grid_comm grid communicator returned by setup
 * @n problem size, there are n sample points along each axis
 * @a value of coefficient a at points x
 * @f value of rhs at points x
 * @v candidate solution at points x
 * @res the residual norm
 */
//void poisson_residual(MPI_Comm &grid_comm, int n, std::vector<real_t> &a, std::vector<real_t> &f,
                      //std::vector<real_t> &v, real_t &res);

void laplacian(const MPI_Comm &grid_comm, int n, const std::vector<real_t> &a,
                   const std::vector<real_t> &u, std::vector<real_t> &lu);
