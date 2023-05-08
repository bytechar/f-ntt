#pragma once

#include <cstdlib>
#include <type_traits>

/**
 * @brief Custom register. Not to be confused with the C/C++ register keyword.
 * @tparam t register buffer type. integral, default = 32 bit signed integer
 */
template <typename t = uint64_t> struct reg {
private:
  static_assert(std::is_integral<t>::value);

public:
  /**
   * @brief Underlying register buffer
   */
  t val;

  // --------- cctor and dtor -----------
  reg() : val(0){};

  explicit reg(t init_val) : val(init_val) {}

  reg(const reg &rhs) : val(rhs.val) {}

  reg(reg &&rhs) noexcept : val(rhs.val) {}

  ~reg() = default;
  // --------- cctor and dtor -----------

  inline reg &operator=(t rhs) {
    val = rhs;
    return *this;
  }

  inline reg &operator=(const reg &rhs) {
    val = rhs.val;
    return *this;
  }

  [[maybe_unused]] inline bool operator[](uint16_t idx) {
    return (val >> idx) & 0x1;
  }

  [[maybe_unused]] inline t operator()() { return val; }

  // --- register modification functions -----------
  // returns reference to itself to allow chaining

  [[maybe_unused]] inline reg &set(uint16_t bit_idx) {
    val |= (1 << bit_idx);
    return *this;
  }

  [[maybe_unused]] inline reg &reset(uint16_t bit_idx) {
    val &= ~(1 << bit_idx);
    return *this;
  }

  [[maybe_unused]] inline reg &flip(uint16_t bit_idx) {
    val ^= (1 << bit_idx);
    return *this;
  }
};

/**
 * @brief overloaded arithmetic scalar multiplication operator
 * @tparam t register buffer type. integral, default = 32 bit signed integer
 * @tparam r return value type. integral, default = 64 bit signed integer
 * @param _1 register 1
 * @param _2 register 2
 * @return integral
 */
template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator*(const reg<t> &_1, const reg<t> &_2) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val * _2.val;
}

/**
 * @brief overloaded arithmetic scalar multiplication operator
 * @tparam t register buffer type. integral, default = 32 bit signed integer
 * @tparam r return value type. integral, default = 64 bit signed integer
 * @param _1 register 1
 * @param _2 integral. default = 32 bit signed integer
 * @return integral
 */
template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator*(const reg<t> &_1, t _2) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val * _2;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator*(t _2, const reg<t> &_1) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val * _2;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator*=(const reg<t> &_1, const reg<t> &_2) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val * _2.val;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator*=(const reg<t> &_1, t _2) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val * _2;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator*=(t _2, const reg<t> &_1) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val * _2;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator+(const reg<t> &_1, const reg<t> &_2) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val + _2.val;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator+(const reg<t> &_1, t _2) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val + _2;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator+(t _2, const reg<t> &_1) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val + _2;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator+=(const reg<t> &_1, const reg<t> &_2) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val + _2.val;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator+=(const reg<t> &_1, t _2) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val + _2;
}

template <typename t = uint64_t, typename r = uint64_t>
constexpr inline r operator+=(t _2, const reg<t> &_1) {
  static_assert(std::is_integral<t>::value && std::is_integral<r>::value);
  return _1.val + _2;
}
