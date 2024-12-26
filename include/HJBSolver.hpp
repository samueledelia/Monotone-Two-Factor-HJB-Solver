#include "PDESolver.hpp"

constexpr uint16_t HJB_T_DIM = 3;

template<std::floating_point Real>
class HJBSolver : public PDESolver<Real, HJB_T_DIM>
{
public:
    
    HJBSolver(const uint32_t N1, const uint32_t N2, const uint32_t N_tau,
                  std::unique_ptr<TwoAssetMinMaxOption<Real>> option, std::pair<Real, Real>& S_max, bool is_sup);
    
    void solve();

protected:

    Eigen::Tensor<Real, HJB_T_DIM> initBoundaryConditions();

    Real getSigma1();

    Real getSigma2();

    const uint32_t N1_, N2_;
    const double dS1_, dS2_;
    const std::pair<Real, Real> S_max_;
    bool is_sup_;
};
