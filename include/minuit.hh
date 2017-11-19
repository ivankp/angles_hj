#ifndef IVANP_MINUIT_HH
#define IVANP_MINUIT_HH

#include "TMinuit.h"
#include <utility>
// #include <boost/type_traits.hpp>

template <typename F>
class minuit final: public TMinuit {
  F f;
  // static constexpr auto arity = boost::function_traits<F>::arity;

  // template <size_t... I>
  // inline auto eval_impl(Double_t* x,std::index_sequence<I...>) {
  //   f(x[I]...);
  // }

public:
  minuit(unsigned npar, const F& f): TMinuit(npar), f(f) { }
  minuit(unsigned npar, F&& f): TMinuit(npar), f(std::move(f)) { }
  minuit(minuit&& r): TMinuit(r.fNpar), f(std::move(r.f)) { }

  inline Int_t Eval(
    Int_t npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag
  ) {
    // fval = eval_impl(par,std::make_index_sequence<N>{});
    fval = f(par);
    return 0;
  }
};

template <typename F>
inline minuit<F> make_minuit(unsigned npar, F&& f) {
  return { npar, std::forward<F>(f) };
}

#endif
