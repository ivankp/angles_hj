#ifndef IVANP_ORDERED_MAP_HH
#define IVANP_ORDERED_MAP_HH

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>

#include "error.hh"

template <typename T>
struct deref_pred : T {
  using T::T;
  template <typename... Args>
  inline auto operator()(const Args&... args) const {
    return T::operator()(*args...);
  }
};

template <typename Iterator>
class deref_iterator {
  Iterator it;
public:
  constexpr deref_iterator(Iterator it): it(it) { }
  inline auto& operator*() noexcept { return **it; }
  inline auto& operator->() noexcept { return *it; }
  inline auto& operator++() noexcept(noexcept(++it)) { ++it; return *this; }
  constexpr bool operator==(const deref_iterator& r) const noexcept
  { return it == r.it; }
  constexpr bool operator!=(const deref_iterator& r) const noexcept
  { return it != r.it; }
  constexpr bool operator<(const deref_iterator& r) const noexcept
  { return it < r.it; }
  constexpr bool operator>(const deref_iterator& r) const noexcept
  { return it > r.it; }
  constexpr Iterator underlying() const noexcept { return it; }
};

template <typename T, typename Key = std::string,
          typename Hash = std::hash<Key>,
          typename KeyEqual = std::equal_to<Key> >
class ordered_map {
  std::unordered_map<Key,T,Hash,KeyEqual> map;
  std::vector<typename decltype(map)::iterator> order;

public:
  using iterator =
    deref_iterator<typename decltype(order)::iterator>;
  using const_iterator =
    deref_iterator<typename decltype(order)::const_iterator>;

  template <typename K>
  auto& operator[](K&& key) {
    auto emp = map.emplace(
      std::piecewise_construct,
      std::forward_as_tuple(std::forward<K>(key)),
      std::forward_as_tuple());
    if (emp.second) order.push_back(emp.first);
    return emp.first->second;
  }
  const T& operator[](const Key& key) const {
    try {
      return map.at(key);
    } catch (const std::out_of_range&) {
      throw ivanp::error("no key \"",key,'\"');
    }
  }
  template <typename... Args1, typename... Args2>
  inline bool emplace(
    std::tuple<Args1...>&& args1, std::tuple<Args2...>&& args2={}
  ) {
    auto emp = map.emplace(
      std::piecewise_construct,
      std::move(args1), std::move(args2));
    if (emp.second) order.push_back(emp.first);
    return emp.second;
  }

  inline iterator begin() noexcept { return order.begin(); }
  inline iterator   end() noexcept { return order.  end(); }
  inline const_iterator begin() const noexcept { return order.begin(); }
  inline const_iterator   end() const noexcept { return order.  end(); }
  inline auto size() const noexcept { return order.size(); }
  inline auto&  back() { return *order. back(); }
  inline auto& front() { return *order.front(); }
  inline const auto&  back() const { return *order. back(); }
  inline const auto& front() const { return *order.front(); }

  void sort() {
    using it_t = typename decltype(order)::value_type;
    std::sort( order.begin(), order.end(),
      [](const it_t& a, const it_t& b) { return a->first < b->first; }
    );
  }
  template <typename Pred>
  void sort(Pred&& pred) {
    using it_t = typename decltype(order)::value_type;
    std::sort( order.begin(), order.end(),
      [&pred](const it_t& a, const it_t& b) { return pred(*a,*b); });
  }

  bool erase_key(const Key& key) {
    const auto it = map.find(key);
    if (it==map.end()) return false;
    order.erase(std::find(order.begin(),order.end(),it));
    map.erase(it);
    return true;
  }
  iterator erase(iterator it) {
    auto u = it.underlying();
    map.erase(*u);
    return order.erase(u);
  }
};

#endif
