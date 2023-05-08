#pragma once
#include "matrix.hpp"
#include "pointwise_matrix_multiply.hpp"
#include "sys_arr.hpp"
#include "utils.hpp"

using sys_arr_type = systolic_array<false, 1024, 1024>;
using mult_array_type1 = mult_array<true, 1024, 1024>;
using mult_array_type2 = mult_array<false, 1024, 1024>;

/**
 *@brief Perform the NTT given the input vector
 *
 * @tparam input_element_type
 * @tparam transform_element_type
 * @param input_vector
 * @return transformation<transform_element_type>
 */
template <typename input_element_type = uint64_t,
          typename transform_element_type = uint64_t>
transformation<transform_element_type>
ntt(const std::vector<input_element_type> &input_vector) {
  transformation<transform_element_type> t;
  // step 0: get n as the size of the vector
  const uint64_t n = input_vector.size();
  // step 1: calculate the minimum working modulus s.t. all numbers belong to
  // [0,M) and 1<=n<M
  const uint64_t max_element =
      *std::max_element(input_vector.begin(), input_vector.end());
  const uint64_t M = max_element >= n ? (max_element + 1) : (n + 1);

#if defined(DEBUG_NTT)
  printf("n: %llu, M: %llu\n", n, M);
#endif
  // step 2: choose k, where k>=1, and N=kn+1, where N = working modulus. N>=M
  t.omega.N = M, t.omega.k = 1;
  for (;; ++t.omega.k) {
    t.omega.N = (t.omega.k * n) + 1;
    if (t.omega.N >= M && is_prime(t.omega.N))
      break;
  }
#if defined(DEBUG_NTT)
  printf("N: %llu, k: %llu\n", t.omega.N, t.omega.k);
#endif
  // step 3: N is prime => multiplicative group of Z has size phi(N) = kn = N-1
  //         Also, there is at least 1 generator g, which is also a primitive
  //         (N-1)th root of unity
  const auto factors_of_N = factorise(t.omega.N - 1);
#if defined(DEBUG_NTT)
  for (auto &factor : factors_of_N)
    printf("Factors of N: %llu ", factor.first);
  printf("\n");
#endif

  uint64_t generator = 0;
  for (decltype(t.omega.N) a = 2; a < t.omega.N; ++a) {
    for (auto &factor : factors_of_N) {
      const uint64_t pow_a = (t.omega.N - 1) / factor.first;
      const uint64_t mod_exp_a = mod_exp(a, pow_a, t.omega.N);
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

  // primitive n-th root of unity, in mod form
  t.omega.W = mod_exp(generator, t.omega.k, t.omega.N);
#if defined(DEBUG_NTT)
  printf("Generator: %llu, Omega: %llu\n", generator, t.omega.W);
#endif

  // generate a nxn unit matrix
  matrix<uint64_t> omega_mtx(n, std::vector<uint64_t>(n, 1));

  for (uint16_t row = 1; row < n; ++row)
    for (uint16_t col = 1; col < n; ++col)
      omega_mtx[row][col] = mod_exp(t.omega.W, uint64_t(row * col), t.omega.N);

#if defined(DEBUG_NTT)
  /*
    std::cout << "-------------------\n";
    for (auto &row : omega_mtx) {
      for (auto &elem : row) {
        std::cout << elem << " ";
      }
      std::cout << "\n";
    }

  */
  std::cout << "--------Input Vector--------\n";
  for (auto &elem : input_vector)
    std::cout << elem << " ";
  std::cout << "\n";
#endif
  sys_arr_type array;
  t.transform = multiply(array, omega_mtx, input_vector
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                         ,
                         t.omega.N
#endif
  );
#if defined(DEBUG_NTT)
  std::cout << "-------Transform NTT_1-D------------\n";
  for (auto &elem : t.transform)
    std::cout << elem << " ";
  std::cout << "-------------------\n";
#endif
  return t;
}

template <typename input_element_type = uint64_t,
          typename transform_element_type = uint64_t>
std::vector<transform_element_type>
inv_ntt(const transformation<input_element_type> &transform) {
  uint64_t omega_inv = 0;
  uint64_t n_inv = 0;

  for (uint64_t a = 1; a < transform.omega.N; ++a) {
    auto a_mod = (a % transform.omega.N);
    if (((a_mod * transform.omega.W) % transform.omega.N) == 1) {
      omega_inv = a;
      break;
    }
  }

  for (uint64_t a = 1; a < transform.omega.N; ++a) {
    auto a_mod = (a % transform.omega.N);
    if ((a_mod * (transform.transform.size() % transform.omega.N)) %
            transform.omega.N ==
        1) {
      n_inv = a;
      break;
    }
  }

#if defined(DEBUG_NTT)
  printf("Omega inverse: %llu, n inverse: %llu\n", omega_inv, n_inv);
#endif

  matrix<uint64_t> omega_mtx(
      transform.transform.size(),
      std::vector<uint64_t>(transform.transform.size(), 1));
  for (uint16_t row = 1; row < transform.transform.size(); ++row)
    for (uint16_t col = 1; col < transform.transform.size(); ++col)
      omega_mtx[row][col] =
          mod_exp(omega_inv, uint64_t(row * col), transform.omega.N);
  /*
  #if defined(DEBUG_NTT)
    std::cout << "-------------------\n";
    for (auto &row : omega_mtx) {
      for (auto &elem : row) {
        std::cout << elem << " ";
      }
      std::cout << "\n";
    }
  #endif
  */

  sys_arr_type array;
  auto result = multiply(array, omega_mtx, transform.transform
#if defined(USE_MAC_VERSION_3)
                         ,
                         transform.omega.N
#endif
  );
  for (auto &elem : result)
    elem = (elem * n_inv) % transform.omega.N;
#if defined(DEBUG_NTT)
  std::cout << "--------INV_Transform---------\n";
  for (auto &elem : result) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;
#endif
  return result;
}

template <typename vector_element_type = uint64_t,
          typename matrix_element_type = uint64_t, typename basetype = uint16_t,
          typename rowtype = uint16_t, typename coltype = uint16_t,
          typename modtype = uint64_t>
std::vector<vector_element_type>
col_ntt(matrix<matrix_element_type> CM1_mtx,
        matrix<matrix_element_type> CM2_mtx, matrix<matrix_element_type> Wm,
        std::vector<vector_element_type> Xbci, basetype b, rowtype N3,
        modtype mod, sys_arr_type &sys_arr, mult_array_type2 &mult_arr) {

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Xbci_vec-----------\n";
  for (auto &elem : Xbci) {
    std::cout << elem << " ";
  }
  std::cout << "\n";
#endif

  // Factorize the column vector into bxN1 matrix
  matrix<matrix_element_type> Xbci_mtx =
      vector_to_matrix<vector_element_type, uint16_t, uint16_t>(Xbci, b,
                                                                N3 / b);

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Xbci_mtx-----------\n";
  for (auto &row : Xbci_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------CM1_mtx-----------\n";
  for (auto &row : CM1_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // Perfrom systolic matrix multiplication
  auto CM1_Xbci_mtx = multiply(sys_arr, CM1_mtx, Xbci_mtx
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                               ,
                               mod
#endif
  );

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------CM1_Xbci_mtx-----------\n";
  for (auto &row : CM1_Xbci_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Wm-----------\n";
  for (auto &row : Wm) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // Perform Pointwise multiplication with Wm array
  auto Ybci_mtx = coeff_multiply(mult_arr, CM1_Xbci_mtx, Wm
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                                 ,
                                 mod
#endif
  );

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Ybci_mtx-----------\n";
  for (auto &row : Ybci_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // Perfrom systolic multiplication with Cm2 matrix
  auto Zbci_mtx = multiply(sys_arr, CM2_mtx, transpose(Ybci_mtx)
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                                                 ,
                           mod
#endif
  );

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Zbci_mtx-----------\n";
  for (auto &row : Zbci_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif
  auto Zbci = matrix_to_vector(Zbci_mtx);

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Zbci-----------\n";
  for (auto &elem : Zbci) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;
#endif

  return Zbci;
};

template <typename vector_element_type = uint64_t,
          typename matrix_element_type = uint64_t, typename basetype = uint16_t,
          typename rowtype = uint16_t, typename coltype = uint16_t,
          typename modtype = uint64_t>
std::vector<vector_element_type>
row_ntt(matrix<matrix_element_type> CM1_mtx,
        matrix<matrix_element_type> CM2_mtx, matrix<matrix_element_type> Wm,
        std::vector<vector_element_type> Xbri, basetype b, rowtype N4,
        modtype mod, sys_arr_type &sys_arr, mult_array_type2 &mult_arr) {

  // Factorize the column vector into N1 row and N2 col
  auto Xbri_mtx = vector_to_matrix<vector_element_type, uint16_t, uint16_t>(
      Xbri, b, N4 / b);

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Xbri_mtx-----------\n";
  for (auto &row : Xbri_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------CM1_mtx-----------\n";
  for (auto &row : CM1_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // Perfrom systolic matrix multiplication
  auto CM1_Xbri_mtx = multiply(sys_arr, CM1_mtx, Xbri_mtx
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                               ,
                               mod
#endif
  );

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------CM1_Xbri_mtx-----------\n";
  for (auto &row : CM1_Xbri_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Wm-----------\n";
  for (auto &row : Wm) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // Perform Pointwise multiplication with Wm array
  auto Ybri_mtx = coeff_multiply(mult_arr, CM1_Xbri_mtx, Wm
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                                 ,
                                 mod
#endif
  );

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Ybri_mtx-----------\n";
  for (auto &row : Ybri_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // Perfrom systolic multiplication with Cm2 matrix
  auto Zbri_mtx = multiply(sys_arr, CM2_mtx, transpose(Ybri_mtx)
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                                                 ,
                           mod
#endif
  );

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Zbri_mtx-----------\n";
  for (auto &row : Zbri_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif
  auto Zbri = matrix_to_vector(Zbri_mtx);

#if defined(DEBUG_NTT_DEEP)
  std::cout << "--------Zbri-----------\n";
  for (auto &elem : Zbri) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;
#endif

  return Zbri;
};

template <typename transform_element_type = uint64_t,
          typename matrix_element_type = uint64_t, typename ntype = uint64_t,
          typename rowtype = uint16_t, typename coltype = uint16_t,
          typename basetype = uint16_t>
transformation<transform_element_type>
ntt_2d(std::vector<transform_element_type> input, rowtype N3, coltype N4,
       basetype b) {

#if defined(DEBUG_NTT)
  std::cout << "----------Input_VEC---------\n";
  for (auto &elem : input) {
    std::cout << elem << " ";
  }
  std::cout << "\n";
#endif
  // 1. get the input vector and factorize it into N3 row and N4 col matrix
  auto input_mtx = vector_to_matrix(input, N3, N4);

#if defined(DEBUG_NTT)
  std::cout << "----------Input_MTX---------\n";
  for (auto &row : input_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif
  // also obtain the matrix transpose
  auto input_mtx_T = transpose(input_mtx);

  // Instantiate the necessary transformation parameters
  transformation<transform_element_type> TN;

  // Instantiate the different primitive roots
  omega_ntt<transform_element_type> WNb;
  omega_ntt<transform_element_type> WN3;
  omega_ntt<transform_element_type> WN4;

  ntype n = input.size();

  // Calculate the necessary transformation parameters (this stores result and
  // WN)
  TN.omega = primitive_omega_ntt_gen(input, n);

  // calculate the primitive roots given the primiary primitive nth root
  WNb = primitive_omega_finder<transform_element_type>(TN.omega, b);
  WN3 = primitive_omega_finder<transform_element_type>(TN.omega, N3);
  WN4 = primitive_omega_finder<transform_element_type>(TN.omega, N4);

  // Compute the WN, Wb, WN3 and WN4 matrix
  auto WN = omega_matrix_gen(TN.omega, N3, N4);
  auto WN3_mtx = omega_matrix_gen(WN3, N3 / b, N3 / b);
  auto WN4_mtx = omega_matrix_gen(WN4, N4 / b, N4 / b);

#if defined(DEBUG_NTT)
  std::cout << "----------WN---------\n";
  for (auto &row : WN) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

#if defined(DEBUG_NTT)
  std::cout << "----------WN3---------\n";
  for (auto &row : WN3_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

#if defined(DEBUG_NTT)
  std::cout << "----------WN4---------\n";
  for (auto &row : WN3_mtx) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // obtain CB matrix from TNb used to calculate CM2 matrix
  auto CB = omega_matrix_gen(WNb, b, b);

#if defined(DEBUG_NTT)
  std::cout << "----------CB---------\n";
  for (auto &row : CB) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // opbtain CB transpose for CM1 calculations
  auto CB_T = transpose(CB);

#if defined(DEBUG_NTT)
  std::cout << "----------CB Transpose---------\n";
  for (auto &row : CB_T) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // Generrate CM1 and CM2 matrix for col ntt
  auto CM1_C = CM1_gen(CB_T, N3 / (b * b));
  auto CM2_C = CM2_gen(CB, N3 / (b * b));

#if defined(DEBUG_NTT)
  std::cout << "----------CM1_C---------\n";
  for (auto &row : CM1_C) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

#if defined(DEBUG_NTT)
  std::cout << "----------CM2_C---------\n";
  for (auto &row : CM2_C) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // Generate CM1 and CM2 matrix for row ntt
  auto CM1_R = CM1_gen(CB_T, N4 / (b * b));
  auto CM2_R = CM2_gen(CB, N4 / (b * b));

#if defined(DEBUG_NTT)
  std::cout << "----------CM1_R---------\n";
  for (auto &row : CM1_R) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

#if defined(DEBUG_NTT)
  std::cout << "----------CM2_R---------\n";
  for (auto &row : CM2_R) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // instantiate col_dft systolic and multi hardware class
  sys_arr_type col_sys_arr;
  mult_array_type2 col_mult_arr;

  // instantiate row_dft systolic and multi hardware class
  sys_arr_type row_sys_arr;
  mult_array_type2 row_mult_arr;

  // instantiate Wn multiplier for twiddle multiplication
  mult_array_type2 output_omega_multiplier;

  // compute col_ntt
  auto Zbc_T = input_mtx_T;
  for (coltype c = 0; c < input_mtx_T.size(); ++c) {
    Zbc_T[c] = col_ntt(CM1_C, CM2_C, WN3_mtx, input_mtx_T[c], b, N3
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                       ,
                       TN.omega.N
#endif
                       ,
                       col_sys_arr, col_mult_arr);
  }

#if defined(DEBUG_NTT)
  std::cout << "----------Zbc_T---------\n";
  for (auto &row : Zbc_T) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  auto Zbc = transpose(Zbc_T);

#if defined(DEBUG_NTT)
  std::cout << "----------Zbc---------\n";
  for (auto &row : Zbc) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  auto ZbcW = coeff_multiply(output_omega_multiplier, Zbc, WN
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                             ,
                             TN.omega.N
#endif
  );

#if defined(DEBUG_NTT)
  std::cout << "----------ZbcW---------\n";
  for (auto &row : ZbcW) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  // compute row_ntt
  auto Zbr = ZbcW;
  for (rowtype r = 0; r < ZbcW.size(); ++r) {
    Zbr[r] = row_ntt(CM1_R, CM2_R, WN4_mtx, ZbcW[r], b, N4
#if defined(USE_MAC_VERSION_3) || defined(USE_MAC_VERSION_4)
                     ,
                     TN.omega.N
#endif
                     ,
                     row_sys_arr, row_mult_arr);
  }

#if defined(DEBUG_NTT)
  std::cout << "----------Zbr---------\n";
  for (auto &row : Zbr) {
    for (auto &elem : row) {
      std::cout << elem << " ";
    }
    std::cout << "\n";
  }
#endif

  TN.transform = matrix_to_vector(Zbr);

#if defined(DEBUG_NTT)
  std::cout << "-------Transform NTT_2-D------------\n";
  for (auto &elem : TN.transform)
    std::cout << elem << " ";
  std::cout << "-------------------\n";
#endif

  return TN;
}