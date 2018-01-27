#include <iostream>
#include <vector>
#include <functional>

typedef double real_t;
typedef void matvec_t (const std::vector<real_t>& /*input */,
                             std::vector<real_t>& /*output*/);


/*
 * @comm MPI communicator (to compute residual)
 * @mv the matvec
 * @diag diagonal component of the matvec
 * @rhs the right hand side
 * @u0 initial guess
 * @eps_r relative tolerance
 * @eps_a absolute tolerance
 * @m check residual every m iterations
 * @k_max the maximum number of iteration
 * @u_final return value
 * @k the number of iteration taken
 */
void jacobi_solver(const MPI_Comm &comm, std::function<matvec_t> mv,
                 std::vector<real_t> &diag, std::vector<real_t> &rhs,
                 std::vector<real_t> &u0, real_t eps_r, real_t eps_a, int k_max,
                 std::vector<real_t> &u_final, int &k);
