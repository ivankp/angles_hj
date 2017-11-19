#ifndef IVANP_TERMCOLOR_MSG_HH
#define IVANP_TERMCOLOR_MSG_HH

#include <iostream>
#include <exception>

#include "termcolor.hpp"
#include "catstr.hh"

namespace ivanp {

template <typename T, typename... TT>
inline void info(const T& x, const TT&... xx) {
  using namespace termcolor;
  std::cout << blue << bold << x << reset << ':';
  if (sizeof...(TT)) {
    std::cout << ' ';
    detail::cat_impl(std::cout,xx...);
  }
  std::cout << std::endl;
}

template <typename T, typename... TT>
inline void error(const T& x, const TT&... xx) {
  using namespace termcolor;
  std::cerr << red << bold << x << reset << ':';
  if (sizeof...(TT)) {
    std::cerr << ' ';
    detail::cat_impl(std::cerr,xx...);
  }
  std::cerr << std::endl;
}

}

std::ostream& operator<<(std::ostream& out, const std::exception& e) {
  using namespace termcolor;
  return out << red << bold << e.what() << reset;
}

#endif
