#include "PDESolver.hpp"

template<std::floating_point Real>
PDESolver<Real>::PDESolver(size_t N1, size_t N2, size_t N_tau, Option<Real>& option)
:option_(option)
{
    U_ = Eigen::Tensor<Real, 3>(N_tau, N1, N2);
}

template<std::floating_point Real>
Eigen::Tensor<Real, 3> PDESolver<Real>::getU()
{
    return U_;
}

template class PDESolver<double>;