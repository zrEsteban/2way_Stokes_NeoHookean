cmake_minimum_required(VERSION 3.16)
project(dealii_pdms_solid LANGUAGES CXX)

find_package(deal.II 9.3 REQUIRED
  HINTS /usr/share/cmake/deal.II /usr/local/share/cmake/deal.II)
deal_ii_initialize_cached_variables()

add_executable(dealiiPdmsSolid dealiiPdmsSolid.cc)
deal_ii_setup_target(dealiiPdmsSolid)

add_executable(testNeoHookean tests/test_material.cc)
target_include_directories(testNeoHookean PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
deal_ii_setup_target(testNeoHookean)

add_executable(testTimeIntegration tests/test_time_integration.cc)
target_include_directories(testTimeIntegration PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(testSignConvention tests/test_sign_convention.cc)

add_executable(testInterfaceSparsity tests/test_interface_sparsity.cc)
target_include_directories(testInterfaceSparsity PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
deal_ii_setup_target(testInterfaceSparsity)

add_executable(testRuntimeNewton tests/test_runtime_newton.cc)
target_include_directories(testRuntimeNewton PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
deal_ii_setup_target(testRuntimeNewton)

enable_testing()
add_test(NAME neo_hookean_material COMMAND testNeoHookean)
add_test(NAME time_integration COMMAND testTimeIntegration)
add_test(NAME sign_convention COMMAND testSignConvention)
