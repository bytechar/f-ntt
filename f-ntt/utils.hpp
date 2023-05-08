#pragma once
#include <bitset>
#include <cmath>
#include <math.h>
#include <random>
#include <set>
#include <unordered_map>
#include <vector>
#include<iostream>

template <typename omega_type> struct omega_ntt {
  uint64_t W;
  uint64_t k;
  uint64_t N;
  uint64_t g;
};

template <typename tranform_element_type = uint64_t> struct transformation {
  std::vector<tranform_element_type> transform;
  omega_ntt<tranform_element_type> omega;
};

/**
num	pow		result	result
4	1101	1	    skip
16	110		4	    apply
256	11		4	    apply
429	1		30	    apply
151	0		445
*/

/**
 * @brief Modular exponent using left to right binary method.
 * @description This is largely similar to square and multiply in the sense that
 * instead of "building up" from square and multiply, we go the other way round
 * and divide and root
 * @param num
 * @param pow
 * @param mod
 * @return
 */
template <typename base_type = uint64_t, typename power_type = uint64_t,
          typename mod_type = uint64_t, typename return_type = uint64_t>
return_type mod_exp(base_type num, power_type pow, mod_type mod) {
  static_assert(std::is_integral<base_type>::value &&
                std::is_unsigned<base_type>::value &&
                std::is_integral<power_type>::value &&
                std::is_unsigned<power_type>::value &&
                std::is_integral<mod_type>::value &&
                std::is_unsigned<mod_type>::value &&
                std::is_integral<return_type>::value &&
                std::is_unsigned<return_type>::value);
  if (mod == 1)
    return 0;          // one divides everything, remainder zero
  uint64_t result = 1; // start with a product of 1
  num = num % mod;     // get the mod of the number
  while (pow > 0) {
    // if it is an odd power, then multiply the result with a single power of
    // num
    if (pow & uint64_t(1))
      result = (result * num) % mod;
    // square the number and mod it
    num = (num * num) % mod;
    pow >>= 1; // right shift the exponent by 1
  }
  return result;
}

template <typename num_type = uint64_t, typename mod_type = uint64_t>
bool mod_eq(num_type lhs, num_type rhs, mod_type mod) {
  if (mod == 0)
    return false;
  if (mod == 1)
    return true;
  if (lhs == 0 || rhs == 0)
    return false;
  if (lhs == rhs)
    return true;
  return (lhs % mod) == (rhs % mod);
}

template <typename t = uint64_t> t bit_reverse(t num) {
  // making sure that only non-negative numbers are bit reversed
  static_assert(std::is_integral<t>::value && std::is_unsigned<t>::value);
  uint64_t bit_length = sizeof(t) * 8;
  t ret_val = 0;
  for (short bit_idx = 0; bit_idx < bit_length; ++bit_idx)
    ret_val |=
        static_cast<t>(((num >> bit_idx) & 1) << (bit_length - bit_idx - 1));
  return ret_val;
}

/**
 * @brief Get the binary representation of a number
 *
 * @tparam t unsigned integral type
 * @param num
 * @return std::string string with binary representation of the number. Assumes
 * little endian representation
 */
template <typename t> std::string get_binary_string(t num) {
  static_assert(std::is_integral<t>::value && std::is_unsigned<t>::value);
  std::bitset<sizeof(t) * 8> b(num);
  return b.to_string();
}
/**
 *@brief Optimized sieve of eratosthenes to check if the number is prime.
 *
 * @tparam t
 * @param num
 * @return true
 * @return false
 */
template <typename t> bool is_prime(t num) {
  static_assert(std::is_integral<t>::value && std::is_unsigned<t>::value);
  if (num < 2)
    return false;

  if (num == 2 || num == 3)
    return true;

  if ((num & 0x1) == 0x0) // if the number >2 and is even, it cannot be prime
    return false;

  const auto sqrt_of_num = std::ceil(std::sqrt(num));

  for (uint64_t i = 3; i <= sqrt_of_num; i += 2) {
    if (num % i == 0)
      return false;
  }
  return true;
}

/**
 *@brief Factorization function finds factors and there frequency
 *
 * @tparam num_type
 * @tparam factors_type
 * @param number
 * @return std::unordered_map<factors_type, factors_type>
 */
template <typename num_type = uint64_t, typename factors_type = uint64_t>
std::unordered_map<factors_type, factors_type> factorise(num_type number) {
  std::unordered_map<factors_type, uint64_t> factors;
  if (number != 0) {
    if (number == 1)
      factors.insert({1, 1});
    else {
      while (number % (num_type)2 == 0) {
        number /= (num_type)2;
        if (factors.find(2) == factors.end())
          factors.insert({2, 0});

        ++(factors[2]);
      }
    }

    const num_type number_sqrt = std::ceil(std::sqrt(number));
    for (num_type divisor = 3; divisor < number_sqrt; divisor += 2) {
      while (number % divisor == 0) {
        if (factors.find(divisor) == factors.end())
          factors.insert({divisor, 0});

        ++(factors[divisor]);
        number /= divisor;
      }
    }
    if (number > 2)
      factors.insert({number, 1});
  }
  return factors;
}

/**
 *@brief Random vector generator
 *
 * @tparam t
 * @param num
 * @return std::vector<t>
 */

template <typename t = uint64_t>
std::vector<t> random_vector_gen(uint32_t num) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<t> dist(std::llround(std::pow(2, 0)),
                                        std::llround(std::pow(2, 7)));

  std::vector<t> random_vector(num); // create a vector of size num

  std::generate(random_vector.begin(), random_vector.end(),
                [&]() { return dist(gen); });
  return random_vector;
}

/**
 *@brief Convert a vector to matrix
 *
 * @tparam T
 * @tparam RowCountType
 * @tparam ColCountType
 * @param vec
 * @param rows
 * @param cols
 * @return std::vector<std::vector<T>>
 */
template <typename transform_element_type = uint64_t,
          typename RowCountType = uint8_t, typename ColCountType = uint8_t>
std::vector<std::vector<transform_element_type>>
vector_to_matrix(const std::vector<transform_element_type> &vec,
                 RowCountType rows, ColCountType cols) {
  static_assert(std::is_integral_v<RowCountType> &&
                    std::is_integral_v<ColCountType>,
                "RowCountType and ColCountType must be integral types");

  const auto expected_size =
      static_cast<RowCountType>(rows) * static_cast<ColCountType>(cols);
  if (vec.size() != expected_size) {
    std::cout <<vec.size() << " "<< (int)rows << " "<<(int)cols << std::endl;
    throw std::invalid_argument("Vector size does not match matrix size");
  }

  std::vector<std::vector<transform_element_type>> mat(
      rows, std::vector<transform_element_type>(cols));
  transform_element_type idx = 0;
  for (RowCountType i = 0; i < rows; ++i) {
    for (ColCountType j = 0; j < cols; ++j) {
      mat[i][j] = vec[idx++];
    }
  }
  return mat;
}

/**
 * @brief Convert a matrix to a 1D vector
 *
 * @tparam T
 * @tparam RowCountType
 * @tparam ColCountType
 * @param mat
 * @return std::vector<T>
 */
template <typename transform_element_type = uint64_t,
          typename RowCountType = uint8_t, typename ColCountType = uint8_t>
std::vector<transform_element_type>
matrix_to_vector(const std::vector<std::vector<transform_element_type>> &mat) {
  static_assert(std::is_integral_v<RowCountType> &&
                    std::is_integral_v<ColCountType>,
                "RowCountType and ColCountType must be integral types");

  const RowCountType rows = static_cast<RowCountType>(mat.size());
  const ColCountType cols = static_cast<ColCountType>(mat[0].size());

  std::vector<transform_element_type> vec(rows * cols);
  transform_element_type idx = 0;
  for (const auto &row : mat) {
    if (row.size() != cols) {
      throw std::invalid_argument("Matrix rows have different sizes");
    }
    for (const auto &val : row) {
      vec[idx++] = val;
    }
  }
  return vec;
}

/**
 *@brief Convert a matrix to a vector of matrix
 *
 * @tparam t
 * @tparam rowCount_type
 * @tparam colCount_type
 * @param vec
 * @param r
 * @param c
 * @return std::vector<std::vector<t>>
 */
template <typename t = uint64_t, typename rowCount_type = uint8_t,
          typename colCount_type = uint8_t>
[[maybe_unused]] std::vector<std::vector<std::vector<t>>>
matrix_to_matrixVector(const std::vector<std::vector<t>> &mat, rowCount_type N1,
                       colCount_type N2) {
  if (mat[0].size() != N1 * N2) {
    throw std::invalid_argument(
        "N1 and N2 Factorization invalid for given matrix");
  }
  std::vector<std::vector<std::vector<t>>> result;
  result.reserve(mat.size());
  for (uint8_t row = 0; row < mat.size(); ++row) {
    std::vector<std::vector<t>> rowMatrix(N1, std::vector<t>(N2, mat[row][0]));
    for (uint8_t col = 0; col < mat[row].size(); ++col) {
      rowMatrix[col / N2][col % N2] = mat[row][col];
    }
    result.push_back(std::move(rowMatrix));
  }
  return result;
}

/**
 *@brief CM1 matrix generation
 *
 * @tparam t
 * @tparam repeatCount_type
 * @tparam r
 * @param CB
 * @param repeatCount
 * @return std::vector<std::vector<r>>
 */
template <typename t = uint64_t, typename repeatCount_t = uint8_t,
          typename r = uint64_t>
std::vector<std::vector<r>> CM1_gen(std::vector<std::vector<t>> CB,
                                    repeatCount_t repeatCount) {
  std::vector<std::vector<r>> CM1;
  uint8_t numRows = CB.size();
  CM1.reserve(numRows *
              repeatCount); // pre-allocate memory for better performance
  for (repeatCount_t i = 0; i < repeatCount; ++i) {
    CM1.insert(CM1.end(), CB.begin(), CB.end());
  }
  return CM1;
}
/**
 *@brief CM2 Matrix generation
 *
 * @tparam t
 * @tparam repeatCount_type
 * @tparam r
 * @param CB
 * @param repeatCount
 * @return std::vector<std::vector<r>>
 */
template <typename t = uint64_t, typename repeatCount_t = uint8_t,
          typename r = uint64_t>
std::vector<std::vector<r>> CM2_gen(std::vector<std::vector<t>> CB,
                                    repeatCount_t repeatCount) {
  std::vector<std::vector<r>> CM2;
  uint8_t numRows = CB.size();
  uint8_t numCols = CB[0].size();
  CM2.resize(numRows);
  for (uint8_t i = 0; i < numRows; ++i) {
    CM2[i].reserve(numCols * repeatCount);
    for (repeatCount_t j = 0; j < repeatCount; ++j) {
      CM2[i].insert(CM2[i].end(), CB[i].begin(), CB[i].end());
    }
  }
  return CM2;
}

/**
 *@brief Matrix transpose operation
 *
 * @tparam t
 * @tparam r
 * @param inputMatrix
 * @param outputMatrix
 */
template <typename T>
std::vector<std::vector<T>>
transpose(const std::vector<std::vector<T>> &matrix) {
  if (matrix.empty()) {
    return {};
  }

  const uint8_t rows = matrix.size();
  const uint8_t cols = matrix[0].size();

  std::vector<std::vector<T>> result(cols, std::vector<T>(rows));

  for (uint8_t i = 0; i < rows; ++i) {
    for (uint8_t j = 0; j < cols; ++j) {
      result[j][i] = matrix[i][j];
    }
  }

  return result;
}

/**
 *@brief  Array Multiplier to generate the matrix-vector point wise product
 *
 * @tparam t
 * @tparam r
 * @param v : vector that gets multiplied into matrix
 * @param m : matrix to be pointWise multiplies
 * @return std::vector<r>
 */
template <typename t = uint64_t, typename r = uint64_t,
          typename mod_type = uint64_t>
std::vector<r> pointwise_multiply(std::vector<r> v,
                                  std::vector<std::vector<r>> m, mod_type N) {
  static_assert(std::is_same<t, r>::value,
                "Type of vector and matrix elements must match.");
  static_assert(v.size() == m.size(),
                "Vector size must match matrix row size.");
  std::vector<std::vector<r>> result(m.size(), std::vector<r>(m[0].size()));
  for (uint8_t row = 0; row < m.size(); ++row) {
    for (uint8_t col = 0; col < m[row].size(); ++col) {
      result[row][col]((v[row] * m[row][col]) % N);
    }
  }
  return result;
}

/**
 *@brief Omega nth root in mod form generation
 *
 * @tparam input_element_type
 * @tparam transform_element_type
 * @tparam n_type
 * @param input_vector
 * @param n
 * @return transformation<transform_element_type>
 */
template <typename input_element_type = uint64_t,
          typename transform_element_type = uint64_t,
          typename ntype = uint64_t>
omega_ntt<transform_element_type>
primitive_omega_ntt_gen(const std::vector<input_element_type> &input_vector,
                        ntype n) {

  omega_ntt<transform_element_type> om;
  // step 0: get n as the size of the vector
  // n (provided as input)
  // step 1: calculate the minimum working modulus s.t. all numbers belong to
  // [0,M) and 1<=n<M
  const uint64_t max_element =
      *std::max_element(input_vector.begin(), input_vector.end());
  const uint64_t M = max_element >= n ? (max_element + 1) : (n + 1);

  // step 2: choose k, where k>=1, and N=kn+1, where N = working modulus. N>=M
  om.N = M, om.k = 1;
  for (;; ++om.k) {
    om.N = (om.k * n) + 1;
    if (om.N >= M && is_prime(om.N))
      break;
  }

  // step 3: N is prime => multiplicative group of Z has size phi(N) = kn = N-1
  //         Also, there is at least 1 generator g, which is also a primitive
  //         (N-1)th root of unity
  const auto &factors_of_N = factorise(om.N - 1);

  uint64_t generator = 0;
  for (decltype(om.N) a = 2; a < om.N; ++a) {
    for (auto &factor : factors_of_N) {
      const uint64_t pow_a = (om.N - 1) / factor.first;
      const uint64_t mod_exp_a = mod_exp(a, pow_a, om.N);

#if defined(DEBUG_NTT)
      printf("Checking a=%llu is a generator using factor=%llu, Power: %llu, "
             "mod exp "
             "of a: %llu\n",
             a, factor.first, pow_a, mod_exp_a);
#endif

      if (mod_exp_a != 1)
        generator = a;
      else {
        generator = 0;
        break;
      }
    }
    if (generator > 0)
      break;
  }

  // set the generator
  om.g = generator;

  // primitive n-th root of unity, in mod form
  om.W = mod_exp(generator, om.k, om.N);

#if defined(DEBUG_NTT)
  printf("Primary primtive root:\nn: %llu, Working M: %llu, k: %llu, "
         "Generator: %llu, W: %llu\n",
         n, om.N, om.k, om.g, om.W);
#endif

  return om;
}

template <typename matrix_element_type = uint64_t,
          typename transformation_element_type = uint64_t,
          typename row_element_type = uint8_t,
          typename col_element_type = uint8_t>
matrix<matrix_element_type>
omega_matrix_gen(omega_ntt<transformation_element_type> &o,
                 row_element_type row, col_element_type col) {
  // generate a nxn unit matrix
  matrix<uint64_t> omega_mtx(row, std::vector<uint64_t>(col, 1));

  // take a primitive root and generate a omega matrix of given row and col
  for (row_element_type r = 1; r < row; ++r)
    for (col_element_type c = 1; c < col; ++c)
      omega_mtx[r][c] = mod_exp(o.W, uint64_t(r * c), o.N);

  return omega_mtx;
}

template <typename transfrom_element_type>
omega_ntt<transfrom_element_type>
primitive_omega_finder(omega_ntt<transfrom_element_type> omega,
                       transfrom_element_type b) {

  omega_ntt<transfrom_element_type> om;

  om.g = omega.g;

  // given input M is the wotking modulous
  om.N = omega.N;

  // choose k, where k>=1, and N=kn+1, where N = working modulus. N>=M
  //  this will work since input n (either N3, b, or N4) is a factor of N
  om.k = (omega.N - 1) / b;

  // primitive n-th root of unity, in mod form
  om.W = mod_exp(om.g, om.k, om.N);

#if defined(DEBUG_NTT)
  printf(
      "For base n: %llu, Working M: %llu, k:%llu, Generator: %llu, "
      "W:%llu\n",
      b, om.N, om.k, om.g, om.W);
#endif

  return om;
}
