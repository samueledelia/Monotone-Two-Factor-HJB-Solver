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

    Real getU_element(size_t i, size_t j) const;

    Real getS1(size_t i) const;

    Real getS2(size_t i) const;

    Real computeGamma(Real rho, uint32_t i, uint32_t j) const;

    //std::vector determineOptimalControl(const Eigen::Tensor<Real, HJB_T_DIM>& W, uint32_t i, uint32_t j); 
    
    const uint32_t N1_, N2_;
    Vector S1_, S2_;
    const double dS1_, dS2_;
    const std::pair<Real, Real> S_max_;
    bool is_sup_;
};
