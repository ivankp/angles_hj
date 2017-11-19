#ifndef IVANP_LORENTZ_VECTOR
#define IVANP_LORENTZ_VECTOR

namespace ivanp {

struct lorentz_vector {
  double t, x, y, z;

  [[ gnu::const ]] constexpr lorentz_vector operator+(const lorentz_vector& p) const noexcept {
    return { t+p.t, x+p.x, y+p.y, z+p.z };
  }
  [[ gnu::const ]] constexpr lorentz_vector operator-(const lorentz_vector& p) const noexcept {
    return { t-p.t, x-p.x, y-p.y, z-p.z };
  }
  [[ gnu::const ]] constexpr double operator*(const lorentz_vector& p) const noexcept {
    return t*p.t - (x*p.x + y*p.y + z*p.z);
  }

  lorentz_vector& operator+=(const lorentz_vector& p) noexcept {
    t += p.t; x += p.x; y += p.y; z += p.z;
    return *this;
  }
  lorentz_vector& operator-=(const lorentz_vector& p) noexcept {
    t -= p.t; x -= p.x; y -= p.y; z -= p.z;
    return *this;
  }
  lorentz_vector& operator*=(double c) noexcept {
    t *= c; x *= c; y *= c; z *= c;
    return *this;
  }
  lorentz_vector& operator/=(double c) noexcept {
    t /= c; x /= c; y /= c; z /= c;
    return *this;
  }
};

[[ gnu::const ]] constexpr lorentz_vector operator*(double c, const lorentz_vector& p) noexcept {
  return { c*p.t, c*p.x, c*p.y, c*p.z };
}
[[ gnu::const ]] constexpr lorentz_vector operator*(const lorentz_vector& p, double c) noexcept {
  return { c*p.t, c*p.x, c*p.y, c*p.z };
}
[[ gnu::const ]] constexpr lorentz_vector operator/(const lorentz_vector& p, double c) noexcept {
  return { c/p.t, c/p.x, c/p.y, c/p.z };
}

std::ostream& operator<<(std::ostream& s, const lorentz_vector& p) {
  return s << '(' << p.t <<','<< p.x <<','<< p.y <<','<< p.z << ')';
}

template <unsigned i> constexpr double get(const lorentz_vector& p) noexcept;
template <> [[ gnu::const ]] constexpr double get<0>(const lorentz_vector& p) noexcept { return p.t; }
template <> [[ gnu::const ]] constexpr double get<1>(const lorentz_vector& p) noexcept { return p.x; }
template <> [[ gnu::const ]] constexpr double get<2>(const lorentz_vector& p) noexcept { return p.y; }
template <> [[ gnu::const ]] constexpr double get<3>(const lorentz_vector& p) noexcept { return p.z; }

} // end namespace ivanp

#endif
