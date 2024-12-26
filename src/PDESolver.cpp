#include "PDESolver.hpp"

template <std::floating_point Real, uint16_t T_DIM>
PDESolver<Real, T_DIM>::PDESolver(Option<Real>& option, const uint32_t N_tau)
    :N_tau_(N_tau), option_(option)
{
    U_ = Eigen::Tensor<Real, T_DIM>();
    U_.setZero();
}

template<std::floating_point Real, uint16_t T_DIM>
Eigen::Tensor<Real, T_DIM>& PDESolver<Real, T_DIM>::getU()
{
    return U_;
}

template class PDESolver<double, 2>;
template class PDESolver<double, 3>;
