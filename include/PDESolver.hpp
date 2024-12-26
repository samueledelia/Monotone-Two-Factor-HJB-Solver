#ifndef HH__PDESOLVERSHPP_HH
#define HH__PDESOLVERSHPP_HH
#include <concepts>
#include <memory>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <Option.hpp>

template<std::floating_point Real, uint16_t T_DIM>
class PDESolver {
public:

    PDESolver(const uint32_t N_tau, OptionBase<Real>& option);

    PDESolver(const uint32_t N_tau, std::unique_ptr<OptionBase<Real>> option_ptr);

    Eigen::Tensor<Real, T_DIM>& getU();

protected:
    const uint32_t N_tau_;
    Eigen::Tensor<Real, T_DIM> U_;
    std::unique_ptr<OptionBase<Real>> option_;
};


#endif //HH__PDESOLVERSHPP_HH
