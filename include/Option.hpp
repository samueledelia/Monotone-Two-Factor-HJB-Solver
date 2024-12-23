#ifndef HH_OPTION_HH
#define HH_OPTION_HH
#include <functional>
#include <concepts>

template<std::floating_point Real>
class Option {
public:
    Option(std::move_only_function<Real(Real, Real)>& payoff,
          const Real r, const std::pair<Real, Real>& sigmas_1,
          const std::pair<Real, Real>& sigmas_2, 
          const std::pair<Real, Real>& rho, const Real expiry);

    Real evaluate(const Real S1, const Real S2);

    Real getExpiry();

    Real getDiscountRate();

    std::pair<Real, Real> getSigmas_1();

    std::pair<Real, Real> getSigmas_2();

    std::pair<Real, Real> getRhos();

private:
    std::move_only_function<Real(Real, Real)> payoff; // strike is in the payoff fn
    Real r;
    std::pair<Real, Real> sigmas_1;
    std::pair<Real, Real> sigmas_2;
    std::pair<Real, Real> rho;
    Real expiry;
};


#endif // HH_OPTION_HH