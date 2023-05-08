#pragma once

#include "sys_arr.hpp"
#include <vector>

template <typename t = uint64_t> using matrix = std::vector<std::vector<t>>;

// naive matrix matrix multiplication
template <bool weight_stationary, uint16_t width, uint16_t height,
          typename mtx1_element_type = uint64_t,
          typename mtx2_element_type = uint64_t,
          typename result_element_type = uint64_t>
matrix<result_element_type>
multiply(systolic_array<weight_stationary, width, height> array,
         const matrix<mtx1_element_type> &m1,
         const matrix<mtx2_element_type> &m2
#ifdef USE_MAC_VERSION_3
         ,
         uint64_t mod = 0
#endif
#ifdef USE_MAC_VERSION_4
         ,
         uint64_t mod = 0
#endif
) {
  // make sure that the matrix dimensions make sense
  assert(!m1.empty() && !m2.empty() && (m1[0].size() == m2.size()));

  // intialize the result matrix dimensions
  matrix<result_element_type> result;
  result.resize(m1.size());
  for (auto &m : result)
    m.resize(m2[0].size());

  /**
   * Idea here is that we take a row of m1 and col of m2 and process that and
   * store that result in a block. The way that we locate that block is simple
   * a[i,k] * b[k,j] + prev_product = processing_elem[i,j]
   * Refer to this for a basic algorithm -
   * http://ecelabs.njit.edu/ece459/lab3.php
   * Essentially we achieve the numeric expansion as shown in Figure 3.1
   */
  uint16_t row_idx = 0, col_idx = 0, traversal_idx = 0;
  while (row_idx < m1.size()) // if we reach the last element, break the loop
  {
    // get the processing element
    auto &p = array(row_idx, col_idx);

    // uncomment for debugging
    // std::cout << "Row: " << row_idx << " Col: " << col_idx << " Traversal: "
    // << traversal_idx << "\n"; std::cout << "Before pe call: " << p.rf() <<
    // "\n";

    // apply inputs and store the result
    p.apply(m1[row_idx][traversal_idx], m2[traversal_idx][col_idx]
#ifdef USE_MAC_VERSION_3
            ,
            mod
#endif
#ifdef USE_MAC_VERSION_4
            ,
            mod
#endif
    );

    // collect the result in the correct index of the result matrix
    // ideally in the case of a true Systolic array, we won't need to do this
    // result[row_idx][col_idx] = p.rf();
    result[row_idx][col_idx] = p.template rf<pe_reg::stationary>();

    // uncomment for debugging
    // std::cout << "After pe call: " << p.rf() << "\n";

    /**
     * Idea here is that we take a row of m1 and col of m2 and process that and
     * store that result in a block. The way that we locate that block is simple
     * a[i,k] * b[k,j] + prev_product = processing_elem[i,j]
     * Refer to this for a basic algorithm -
     * http://ecelabs.njit.edu/ece459/lab3.php
     * Essentially we achieve the numeric expansion as shown in Figure 3.1
     */
    if (traversal_idx == (m1[0].size() - 1)) {
      // increment col_idx by 1 and keep the row idx the same
      if (col_idx == (m2[0].size() - 1)) {
        ++row_idx;
        col_idx = 0;
      } else {
        ++col_idx;
      }
      traversal_idx = 0;
      continue;
    } else
      ++traversal_idx;
  }
  return result;
}

// naive matrix-vector assuming m2 is a column vector
template <bool weight_stationary, uint16_t width, uint16_t height,
          typename mtx1_element_type = uint64_t,
          typename mtx2_element_type = uint64_t,
          typename result_element_type = uint64_t>
std::vector<result_element_type>
multiply(systolic_array<weight_stationary, width, height> array,
         const matrix<mtx1_element_type> &m1,
         const std::vector<mtx2_element_type> &m2
#ifdef USE_MAC_VERSION_3
         ,
         uint64_t mod
#endif
#ifdef USE_MAC_VERSION_4
         ,
         uint64_t mod
#endif
) {
  assert(!m1.empty() && !m2.empty() && (m1.size() == m2.size()));
  std::vector<result_element_type> result(m1.size(), 0);
  for (uint64_t row_idx = 0, col_idx = 0; row_idx < m1.size();) {
    // if we reach the last element, break the loop
    if (row_idx >= m1.size())
      break;

    // get the processing element
    auto &p = array(row_idx, 0);

    // apply inputs and store the result
    p.apply(m1[row_idx][col_idx], m2[col_idx]
#ifdef USE_MAC_VERSION_3
            ,
            mod
#endif
#ifdef USE_MAC_VERSION_4
            ,
            mod
#endif
    );
    result[row_idx] = p.template rf<pe_reg::stationary>();

    // collect the result in the correct index of the result matrix
    // ideally in the case of a true Systolic array, we won't need to do this
    // result[row_idx][col_idx] = p.rf();

    /**
     * Idea here is that we take a row of m1 and col of m2 and process that and
     * store that result in a block. The way that we locate that block is simple
     * a[i,k] * b[k,j] + prev_product = processing_elem[i,j]
     * Refer to this for a basic algorithm -
     * http://ecelabs.njit.edu/ece459/lab3.php
     * Essentially we achieve the numeric expansion as shown in Figure 3.1
     */
    if (col_idx == (m2.size() - 1)) {
      ++row_idx;
      col_idx = 0;
    } else
      ++col_idx;
  }
  return result;
}

// naive vector-vector multiplication assuming m1 is a row vector and m2 is a
// column vector
template <bool weight_stationary, uint16_t width, uint16_t height,
          typename mtx1_element_type = uint64_t,
          typename mtx2_element_type = uint64_t,
          typename result_element_type = uint64_t>
result_element_type
multiply(systolic_array<weight_stationary, width, height> array,
         const std::vector<mtx1_element_type> &m1,
         const std::vector<mtx2_element_type> &m2
#ifdef USE_MAC_VERSION_3
         ,
         uint64_t mod
#endif
#ifdef USE_MAC_VERSION_4
         ,
         uint64_t mod
#endif
) {
  assert(!m1.empty() && !m2.empty() && (m1.size() == m2.size()));
  result_element_type result;
  for (uint64_t idx = 0; idx < m1.size(); ++idx) {
    // get the processing element
    auto &p = array(0, 0);

    // apply inputs and store the result
    p.apply(m1[idx], m2[idx]
#ifdef USE_MAC_VERSION_3
            ,
            mod
#endif
#ifdef USE_MAC_VERSION_4
            ,
            mod
#endif
    );

    // collect the result in the correct index of the result matrix
    // ideally in the case of a true Systolic array, we won't need to do this
    // result[row_idx][col_idx] = p.rf();
    result = p.template rf<pe_reg::stationary>();

    /**
     * Idea here is that we take a row of m1 and col of m2 and process that and
     * store that result in a block. The way that we locate that block is simple
     * a[i,k] * b[k,j] + prev_product = processing_elem[i,j]
     * Refer to this for a basic algorithm -
     * http://ecelabs.njit.edu/ece459/lab3.php Essentially we achieve the
     * numeric expansion as shown in Figure 3.1
     */
  }
  return result;
}