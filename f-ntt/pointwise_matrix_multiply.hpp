#pragma once

#include "multiplier_array.hpp"
#include <vector>

// naive matrix-vector pointwise multiplication
template <bool weight_stationary, uint16_t height, uint16_t width,
          typename mtx_element_type = uint64_t,
          typename multiplier_element_type = uint64_t,
          typename result_element_type = uint64_t>
matrix<result_element_type>
coeff_multiply(mult_array<weight_stationary, height, width> array,
               const matrix<mtx_element_type> &m1,
               const matrix<mtx_element_type> &m2
#ifdef USE_MAC_VERSION_3
               ,
               uint64_t mod = 0
#endif
#ifdef USE_MAC_VERSION_4
               ,
               uint64_t mod = 0
#endif
) {
  assert(!m1.empty() && !m2.empty() && (m1.size() == m2.size()) &&
         (m1[0].size() == m2[0].size()));

  // intialize the result matrix dimensions
  matrix<result_element_type> result;
  result.resize(m1.size());
  for (auto &m : result)
    m.resize(m1[0].size());

  for (uint16_t row_idx = 0; row_idx < m1.size(); ++row_idx) {
    for (uint16_t col_idx = 0; col_idx < m1[row_idx].size(); ++col_idx) {
      // get the processing element
      auto &p = array(row_idx, col_idx);

      // apply inputs and store the result

      if (weight_stationary) {
        p.apply(m1[row_idx][col_idx], 0
#ifdef USE_MAC_VERSION_3
                ,
                mod
#endif
#ifdef USE_MAC_VERSION_4
                ,
                mod
#endif
        );
        result[row_idx][col_idx] = p.template rf<pe_reg::forwarding_1>();
      } else {
        p.apply(m1[row_idx][col_idx], m2[row_idx][col_idx]
#ifdef USE_MAC_VERSION_3
                ,
                mod
#endif
#ifdef USE_MAC_VERSION_4
                ,
                mod
#endif
        );
        result[row_idx][col_idx] = p.template rf<pe_reg::stationary>();
      }
    }
  }
  return result;
}