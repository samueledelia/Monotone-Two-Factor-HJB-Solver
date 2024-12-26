#ifndef HH__PDESOLVERSHPP_HH
#define HH__PDESOLVERSHPP_HH
#include <vector>
#include <functional>
#include <concepts>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <Option.hpp>

// ! Numerical HJB PDE Solver Class
template<std::floating_point Real, uint16_t T_DIM>
class PDESolver {
public:

    PDESolver(Option<Real>& option, const uint32_t N_tau);

    Eigen::Tensor<Real, T_DIM>& getU();
    
protected:
    const uint32_t N_tau_;
    Eigen::Tensor<Real, T_DIM> U_;
    Option<Real> option_;
};


#endif //HH__PDESOLVERSHPP_HH
