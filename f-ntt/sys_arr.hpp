#pragma once

#include "processing_element.hpp"

/**
 * @brief Systolic array skeleton
 * @tparam weight_stationary: Defines if the array is weight stationary
 * @tparam width: Defines the width of the array (no. of columns)
 * @tparam height Define the height of the array (no. of rows)
 */
template <bool weight_stationary, uint16_t width, uint16_t height>
class systolic_array {
public:
  // typedef a processing element
  typedef pe<weight_stationary, uint64_t, uint64_t, pe_reg::max>
      pe; // {bool, t, r, pe_reg_type}
private:
  // processing elements. These are stored as a 1-d vector, but can store 2-d
  // structures.
  std::vector<pe> array;
  /**
   * For example a 2-D structure can be stored as a 1-d vector
   *     0  1   2
   * 0 [pe][pe][pe]
   * 1 [pe][pe][pe]
   * 2 [pe][pe][pe]
   *
   *stored as
   *[pe][pe][pe].[pe][pe][pe]..[pe][pe][pe]
   */
public:
  [[maybe_unused]] const uint64_t sys_arr_height = height;
  [[maybe_unused]] const uint64_t sys_arr_width = width;

  systolic_array() noexcept : array(width * height, pe()) {}

  systolic_array(const systolic_array &rhs) : array(rhs.array) {}

  systolic_array(systolic_array &&rhs) noexcept : array(std::move(rhs.array)) {}

  void reset_macs() {
    for (uint16_t i = 0; i < height; ++i) {
      for (uint16_t j = 0; j < width; ++j) {
        array[(width * i) + j].template rf<pe_reg::stationary>(0);
        array[(width * i) + j].template rf<pe_reg::forwarding_1>(0);
        array[(width * i) + j].template rf<pe_reg::forwarding_2>(0);
      }
    }
  }

  /**
   *@brief Overloaed operator gets the pe given at index
   *
   * @param elem
   * @return pe&
   */
  // assuming 1-d structure
  pe &operator()(uint16_t elem) { return array[elem]; }

  // assuming 2-d structure
  pe &operator()(uint16_t row, uint16_t col) {
    return array[(width * row) + col];
  }
};
