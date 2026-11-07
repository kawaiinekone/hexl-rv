cmake_minimum_required(VERSION 3.16)
project(hexl_rv VERSION 0.1.0 LANGUAGES CXX C ASM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(-Wall -Wextra -O3 -fPIC)
include_directories(include)

set(HEXL_RV_SOURCES
    src/eltwise_mult_mod.cpp
    src/eltwise_add_sub_mod.cpp
    src/ntt.cpp
)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "riscv64")
    message(STATUS "Configuring with native RV64 assembly acceleration")
    list(APPEND HEXL_RV_SOURCES src/kernel_rv64/eltwise_mult_mod_rv64.s)
    add_compile_definitions(HEXL_RV_HAS_RV64_ASM=1)
endif()

add_library(hexl_rv STATIC ${HEXL_RV_SOURCES})

enable_testing()

add_executable(test_eltwise_mult_mod tests/test_eltwise_mult_mod.cpp)
target_link_libraries(test_eltwise_mult_mod PRIVATE hexl_rv)
add_test(NAME test_eltwise_mult_mod COMMAND test_eltwise_mult_mod)

add_executable(test_eltwise_add_sub_mod tests/test_eltwise_add_sub_mod.cpp)
target_link_libraries(test_eltwise_add_sub_mod PRIVATE hexl_rv)
add_test(NAME test_eltwise_add_sub_mod COMMAND test_eltwise_add_sub_mod)

add_executable(test_ntt tests/test_ntt.cpp)
target_link_libraries(test_ntt PRIVATE hexl_rv)
add_test(NAME test_ntt COMMAND test_ntt)
