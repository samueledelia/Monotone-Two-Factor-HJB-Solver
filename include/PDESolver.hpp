#ifndef HH__PDESOLVERSHPP_HH
#define HH__PDESOLVERSHPP_HH
#include <concepts>
#include <memory>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <Option.hpp>

template<std::floating_point Real, uint16_t T_DIM>
class PDESolver {
public:

    PDESolver(OptionBase<Real>& option, const uint32_t N_tau);

    PDESolver(std::unique_ptr<OptionBase<Real>> option_ptr, const uint32_t N_tau);

    virtual Eigen::Tensor<Real, T_DIM>& getU();

protected:
    const uint32_t N_tau_;
    Eigen::Tensor<Real, T_DIM> U_;
    std::unique_ptr<OptionBase<Real>> option_;
};


#endif //HH__PDESOLVERSHPP_HH

