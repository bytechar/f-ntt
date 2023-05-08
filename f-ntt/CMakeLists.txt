cmake_minimum_required(VERSION 3.21)
project(ntt_implementation)

set(CMAKE_CXX_STANDARD 17)

add_executable(ntt_implementation
        main.cpp
        utils.hpp
        register.hpp
        processing_element.hpp
        sys_arr.hpp
        matrix.hpp
        )
