#include "poisson.h"
#include "jacobi.h"
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include <functional>
#include <sys/time.h>
#include <stdlib.h>

#define pi 3.14


void jacobi_solver(const MPI_Comm &comm, std::function<matvec_t> mv,
                 std::vector<real_t> &diag, std::vector<real_t> &rhs,
                 std::vector<real_t> &u0, real_t eps_r, real_t eps_a, int k_max,
                 std::vector<real_t> &u_final, int &k)
{
	/* STEPS
	 * 
	 * Call laplacian
	 * Get the residual somehow. Individual to each process
	 * Perform MPI_Reduce to obtain the total residual
	 * Check for the residual value if it falls in the range
	 * If it does then stop, if it does not then iterate again 
	 * to go back to the beginning of the loop and call laplacian
	 *
	 */	

	int r;
#if 0	
	for (int k=0; k<u0.size(); k++)
	{
		printf("%lf\n", u0[k]);
	}
#endif
	
	real_t loc_residual = 0;
	std::vector<real_t> point_residual;
	point_residual.resize(u0.size());
	real_t total_residual = 0;
	real_t sqrt_total_residual = -1;
	real_t initial_residual = 0;

	int ctr = 4;
#if 1
	while(ctr)
	{
		MPI_Comm_rank (comm, &r);
		
		/*
		 * calling laplacian here (matvec function)
		 */
		mv (u0, u_final);
		
		/*
		 * calculate residual here for each processor
		 */	
		
		for (int i=0; i<u_final.size(); i++)
		{
			point_residual[i] = rhs[i] - u_final[i];
			loc_residual += (point_residual[i] * point_residual[i]);
		}

		/*
		 * Waiting for all the processes here
		 * to complete executing.
		 */

		MPI_Barrier (comm);
	
		printf("LOC RESIDUAL %lf %d\n", loc_residual, r);	
		/*
		 * Reduce all the local residual to the root node which
		 * is the 0th processor.
		 */
		MPI_Reduce (&loc_residual, &total_residual, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
		MPI_Barrier (comm);

#if 0
		if (r==0)
			printf("TOT RESIDUAL %lf %d\n", total_residual, r);
#endif

		k++;

#if 1	
		if (r == 0)
		{
			//printf("SQRT TOT RESIDUAL %lf", sqrt_total_residual);
#if 0
			if (sqrt_total_residual == -1)
			{
				sqrt_total_residual = sqrt(total_residual);
				initial_residual = sqrt_total_residual;
				printf("I CAME HERE INITIALLY\n");
			}
			else
#endif
#if 1
			{
				/*
				 * calculate here if the sqrt of the total residual
				 * is within the limits of the expected value.
				 */
				
				sqrt_total_residual = sqrt (total_residual);					
				printf("SQRT TOT RESIDUAL %lf", sqrt_total_residual);
				if ((sqrt_total_residual < ((eps_r * initial_residual) + eps_a)) || (k > k_max))
				{
					//break;
					//MPI_Finalize();
				}
			}
#endif
		}

		/*
		 * Obtaining u+ (the updated next value of u) here
		 */
#if 0
		for (int j=0; j<u0.size(); j++)
		{
			u0[j] = (point_residual[j] + (diag[j]*u_final[j]))/(diag[j] * 1.0) ;	
		}

		MPI_Barrier (comm);
		ctr--;
#endif
#endif
		for (int j=0; j<u0.size(); j++)
		{
			u0[j] = (point_residual[j] + (diag[j]*u0[j]))/(diag[j] * 1.0) ;
			u_final[j] = 0;	
		}
	
		loc_residual = 0;
		total_residual = 0;
		MPI_Barrier (comm);

		ctr--;


	}
#endif

}

int main(int argc, char* argv[])
{
	std::vector<point_t> x;
	MPI_Init (&argc, &argv);

	int r;
	MPI_Comm_rank (MPI_COMM_WORLD, &r);
	
	if (argc < 2) {
		if (r==0) 
		{
			printf ("Usage: %s n\n", argv[0]);
			MPI_Finalize();
			exit(1);
		}
	}

	int n = atoi (argv[1]);

	MPI_Comm grid_comm;
	MPI_Comm start_comm = MPI_COMM_WORLD;
	poisson_setup (start_comm, n, grid_comm, x);

	std::vector<real_t> a(0);
	std::vector<real_t> u(0);
	std::vector<real_t> u_final(0);
	std::vector<real_t> lu(0);
	std::vector<real_t> diag(0);
	std::vector<real_t> f(0);

	real_t eps_r = std::pow (10, -5);
	real_t eps_a = std::pow (10, -12);
	int k_max = 1000;
	int k;

	for (auto p : x)
	{
		a.push_back(12);
		u.push_back(sin(4*pi*p.x) * sin(10*pi*p.y) * sin(14*pi*p.z));
		//printf("%lf %lf %lf %lf\n", p.x, p.y, p.z, sin(4*pi*p.x) * sin(10*pi*p.y) * sin(14*pi*p.z));
		//printf("%0.9lf\n", (real_t)sin(4*pi*p.x) * sin(10*pi*p.y) * sin(14*pi*p.z));
		//diag.push_back(12 - (6 * n * n * 1.0));
		diag.push_back(6.0);
		f.push_back((real_t)sin(2*pi*p.x) * sin(12*pi*p.y) * sin(8*pi*p.z));
	}

	lu.resize (u.size());
	u_final.resize(u.size());

	auto mv = std::bind(laplacian, grid_comm, n, a, std::placeholders::_1, std::placeholders::_2);


	jacobi_solver (grid_comm, mv, 
					diag, f, 
					u, eps_r, eps_a, k_max, 
					u_final, k);	


	MPI_Finalize();

}
