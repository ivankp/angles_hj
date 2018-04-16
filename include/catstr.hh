// Developed by Ivan Pogrebnyak, MSU

#ifndef IVANP_CATSTR_HH
#define IVANP_CATSTR_HH

#include <string>
#include <sstream>

namespace ivanp {

template <typename... T>
inline std::string cat(T&&... x) {
  std::stringstream ss;
  using expander = int[];
  (void)expander{0, ((void)(ss << std::forward<T>(x)), 0)...};
  return ss.str();
}
inline std::string cat() { return { }; }

inline std::string cat(std::string x) { return x; }
inline std::string cat(const char* x) { return x; }

}

#endif
