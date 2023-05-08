#pragma once

#include "register.hpp"
#include <vector>
template <typename t = uint64_t> using reg_file = std::vector<reg<t>>;
using pe_reg_type = uint16_t;

/**
 *@brief Each PE will contain 3 reisters,
 * stationary: if(weight_stationary) then stores the weights, otherwise outputs
 *stired forwarding_1: forwarding_2: sends the outputs to further PE (incase of
 *)
 */
enum class pe_reg : pe_reg_type {
  min = 0,
  stationary = min,
  forwarding_1 = 1,
  forwarding_2 = 2,
  max = forwarding_2 + 1
};

template <bool weight_stationary = true, typename t = uint64_t,
          typename r = uint64_t, pe_reg pe_reg_file_depth = pe_reg::max>
class pe {
  /**
   * @brief Create a SPE register file
   */
  reg_file<r> pe_reg_file;

public:
  /**
   * @brief Default cctor initializes the PE register file with 3 registers for
   * 1-D or a 2-D configuration.
   */
  pe()
      : pe_reg_file(static_cast<pe_reg_type>(
                        pe_reg_file_depth), // size of register file
                    reg(static_cast<r>(weight_stationary ? 1 : 0))) // reg type,
  {}

  /**
   * @brief Copy cctor
   * @param rhs to copy from
   */
  pe(const pe &rhs) : pe_reg_file(rhs.pe_reg_file) {}

  /**
   * @brief Move cctor
   * @param rhs to move from
   */
  pe(pe &&rhs) noexcept : pe_reg_file(std::move(rhs.pe_reg_file)) {}

  /**
   * @brief Set data in a particular register on register file
   * @tparam i idx of the register
   * @param val value
   * @return Reference to itself to allow chaining
   */
  template <pe_reg i = pe_reg::stationary> inline pe &rf(r val) {
    pe_reg_file[static_cast<pe_reg_type>(i)] = val;
    return *this;
  }

  /**
   * @brief Get the data stored
   * @tparam i
   * @return
   */
  template <pe_reg i = pe_reg::stationary> [[nodiscard]] inline r rf() const {
    return pe_reg_file[static_cast<pe_reg_type>(i)].val;
  }

#if defined(USE_MAC_VERSION_1)
  /**
   * @brief The "ALU"
   * @param a input 1
   * @param b input 2
   */
  inline void apply(t a, t b) {
    if (weight_stationary) {
      // for weight stationary, we hold
      // 1. weight in r[mac_reg_idx::stationary]
      // 2. bypass input in r[mac_reg_idx::forwarding_2]
      // 3. convolved result in r[mac_reg_idx::forwarding_1]
      rf<pe_reg::forwarding_1>((rf<pe_reg::stationary>() * a) + b);
      rf<pe_reg::forwarding_2>(a);
    } else // output stationary
    {
      // for output stationary, we hold
      // 1. convolved result in r[mac_reg_idx::stationary]
      // 2. bypass input 1 in r[mac_reg_idx::forwarding_1]
      // 3. bypass input 2 in r[mac_reg_idx::forwarding_2]
      rf<pe_reg::stationary>((a * b) + rf<pe_reg::stationary>());
      rf<pe_reg::forwarding_1>(a);
      rf<pe_reg::forwarding_2>(b);
    }
  }
#elif defined(USE_MAC_VERSION_2)
  inline void apply(t a, t b) {
    if (weight_stationary) {
      // output propagating
      rf<pe_reg::forwarding_1>(a * b);
      rf<pe_reg::forwarding_2>(a);
    } else {
      // output stationary
      rf<pe_reg::stationary>((a * b) + rf<pe_reg::stationary>());
      rf<pe_reg::forwarding_1>(a);
      rf<pe_reg::forwarding_2>(b);
    }
  }
#elif defined(USE_MAC_VERSION_3)
  inline void apply(t a, t b, t mod) {
    if (weight_stationary) {
      // output propagating
      rf<pe_reg::forwarding_1>((a * b) % mod);
      rf<pe_reg::forwarding_2>(a);
    } else {
      // output stationary
      rf<pe_reg::stationary>(
          (((a * b) % mod) + (rf<pe_reg::stationary>() % mod)) % mod);
      rf<pe_reg::forwarding_1>(a);
      rf<pe_reg::forwarding_2>(b);
    }
  }

  // when inputs flow from just one direction
#elif defined(USE_MAC_VERSION_4)
  inline void apply(t a, t b, t mod) {
    if (weight_stationary) {
      // output propagating
      // propogate the output to the next pe
      rf<pe_reg::forwarding_1>((((a * rf<pe_reg::stationary>()) % mod) + b) %
        mod);
      rf<pe_reg::forwarding_2>(a);
    }
      else {
        // output stationary
        rf<pe_reg::stationary>(
            (((a * b) % mod) + (rf<pe_reg::stationary>() % mod)) % mod);
        rf<pe_reg::forwarding_1>(a);
        rf<pe_reg::forwarding_2>(b);
      }
    }
#endif
};
