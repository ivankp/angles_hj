#ifndef IVANP_RANDOM_HH
#define IVANP_RANDOM_HH

#include <random>

template <typename F, typename X=std::uniform_real_distribution<double>>
struct function_sampling_distribution {
private:
  F f;
  X _x;
  using x_type = typename X::result_type;
  std::uniform_real_distribution<decltype(f(std::declval<x_type>()))> _y;
  using y_type = typename decltype(_y)::result_type;

public:
  template <typename _F>
  function_sampling_distribution(x_type a, x_type b, y_type max, _F&& f)
  : f(std::forward<_F>(f)), _x(a,b), _y(0,max) { }

  using result_type = x_type;

  template <typename URNG>
  result_type operator()(URNG& g) {
    for (;;) {
      const x_type x = _x(g);
      if (_y(g) < f(x)) return x;
    }
  }
};

template <typename M, typename F>
inline auto sample(double a, double b, M max, F&& f)
-> function_sampling_distribution<F> {
  return { a, b, max, std::forward<F>(f) };
}

#endif
