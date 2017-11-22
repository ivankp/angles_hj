#ifndef IVANP_TERMCOLOR_MSG_HH
#define IVANP_TERMCOLOR_MSG_HH

#include <iostream>
#include <exception>

#include "termcolor.hpp"
#include "catstr.hh"

namespace ivanp {

template <typename Color, typename T, typename... TT>
inline void tc_msg(const Color& color, const T& x, const TT&... xx) {
  std::cout << color << termcolor::bold << x << termcolor::reset << ':';
  if (sizeof...(TT)) {
    std::cout << ' ';
    detail::cat_impl(std::cout,xx...);
  }
  std::cout << std::endl;
}

template <typename... TT>
inline void info(const TT&... xx) { tc_msg(termcolor::blue,xx...); }

template <typename... TT>
inline void warning(const TT&... xx) { tc_msg(termcolor::yellow,xx...); }

template <typename... TT>
inline void error(const TT&... xx) { tc_msg(termcolor::red,xx...); }

}

std::ostream& operator<<(std::ostream& out, const std::exception& e) {
  using namespace termcolor;
  return out << red << bold << e.what() << reset;
}

#endif
