#ifndef HH_OPTION_HH
#define HH_OPTION_HH
#include "OptionBase.hpp"


template<std::floating_point Real>
class OneAssetOption : public OptionBase<Real>
{
public:
    OneAssetOption(std::function<Real(Real)> payoff,
          const Real r, const Real sigma,
          const Real expiry)
    :OptionBase<Real>(r, expiry), payoff_(std::move(payoff)),
    sigma_(sigma) {}

    Real getSigma() override { return sigma_; }

protected:
    Real evaluateOption_(Real S) override {
        return payoff_(S);
    }

private:
    std::function<Real(Real)> payoff_;
    Real sigma_;
};

template class OneAssetOption<double>;

template<std::floating_point Real>
class TwoAssetOption : public OptionBase<Real>
{
public:
    TwoAssetOption(std::function<Real(Real, Real)> payoff,
          const Real r, const Real sigma_1,
          const Real sigma_2, 
          const Real rho, const Real expiry)
    :OptionBase<Real>(r, expiry), payoff_(std::move(payoff)),
    sigma_1_(sigma_1), sigma_2_(sigma_2), rho_(rho) {}

    Real getSigma_1() { return sigma_1_; }

    Real getSigma_2() { return sigma_2_; }

    Real getRho() { return rho_; }

protected:

    Real evaluateOption_(Real S1, Real S2) override {
        return payoff_(S1, S2);
    }

private:
    std::function<Real(Real, Real)> payoff_;
    Real sigma_1_, sigma_2_, rho_;
};

template class TwoAssetOption<double>;

template<std::floating_point Real>
class TwoAssetMinMaxOption : public OptionBase<Real>
{
public:
    TwoAssetMinMaxOption(std::function<Real(Real, Real)> payoff,
          const Real r, const std::pair<Real, Real>& sigmas_1,
          const std::pair<Real, Real>& sigmas_2, 
          const std::pair<Real, Real>& rhos, const Real expiry)
    :OptionBase<Real>(r, expiry), payoff_(std::move(payoff)),
    sigmas_1_(sigmas_1), sigmas_2_(sigmas_2), rhos_(rhos) {}

    Real evaluate(Real S1, Real S2) {return payoff_(S1, S2); }

    std::pair<Real, Real> getSigmas_1() override { return sigmas_1_; }

    std::pair<Real, Real> getSigmas_2() override { return sigmas_2_; }

    std::pair<Real, Real> getRhos() { return rhos_; }
protected:

    Real evaluateOption_(Real S1, Real S2) override {
        return payoff_(S1, S2);
    }

private:
    std::function<Real(Real, Real)> payoff_;
    std::pair<Real, Real> sigmas_1_, sigmas_2_, rhos_;
};

template class TwoAssetMinMaxOption<double>;


#endif // HH_OPTION_HH