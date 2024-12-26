#include "PDESolver.hpp"

template <std::floating_point Real, uint16_t T_DIM>
PDESolver<Real, T_DIM>::PDESolver(const uint32_t N_tau, std::unique_ptr<OptionBase<Real>> option_ptr)
    :N_tau_(N_tau), option_(std::move(option_ptr))
{
    U_ = Eigen::Tensor<Real, T_DIM>();
    U_.setZero();
}

template <std::floating_point Real, uint16_t T_DIM>
PDESolver<Real, T_DIM>::PDESolver(const uint32_t N_tau, OptionBase<Real>& option)
    : option_(N_tau_(N_tau), std::make_unique<OptionBase<Real>>(option))
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
