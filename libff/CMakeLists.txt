set(FF_EXTRASRCS)
set(FF_EXTRALIBS)
if(${CURVE} STREQUAL "BN128")
  set(
    FF_EXTRASRCS

    ${FF_EXTRASRCS}
    algebra/curves/bn128/bn128_g1.cpp
    algebra/curves/bn128/bn128_g2.cpp
    algebra/curves/bn128/bn128_gt.cpp
    algebra/curves/bn128/bn128_init.cpp
    algebra/curves/bn128/bn128_pairing.cpp
    algebra/curves/bn128/bn128_pp.cpp
  )
  set(
    FF_EXTRALIBS

    ${FF_EXTRALIBS}
    zm
  )
endif()

file(
  GLOB_RECURSE
  LIBFF_SOURCE

  algebra/curves/alt_bn128/*.?pp
  algebra/curves/alt_bn128/*.tcc
  algebra/curves/edwards/*.?pp
  algebra/curves/edwards/*.tcc
  algebra/curves/mnt/mnt4/*.?pp
  algebra/curves/mnt/mnt4/*.tcc
  algebra/curves/mnt/mnt46_common.?pp
  algebra/curves/mnt/mnt46_common.tcc
  algebra/curves/mnt/mnt6/*.?pp
  algebra/curves/mnt/mnt6/*.tcc
  algebra/curves/bls12_377/*.?pp
  algebra/curves/bls12_377/*.tcc
  algebra/curves/bw6_761/*.?pp
  algebra/curves/bw6_761/*.tcc
  common/double.cpp
  common/profiling.cpp
  common/utils.cpp
)

add_library(
  ff
  STATIC

  ${LIBFF_SOURCE}
  ${FF_EXTRASRCS}
)
target_link_libraries(
  ff

  GMP::gmp
  ${PROCPS_LIBRARIES}
  ${FF_EXTRALIBS}
)
target_include_directories(
  ff
  PUBLIC ..
)

install(
  DIRECTORY "" DESTINATION "include/libff"
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.tcc"
  PATTERN "tests" EXCLUDE
  PATTERN "examples" EXCLUDE
)

install(
  TARGETS ff DESTINATION lib
)

# Tests
if ("${IS_LIBFF_PARENT}")
  add_executable(
    algebra_bilinearity_test
    EXCLUDE_FROM_ALL

    algebra/curves/tests/test_bilinearity.cpp
  )
  target_link_libraries(
    algebra_bilinearity_test

    ff
    gtest_main
  )

  add_executable(
    algebra_groups_test
    EXCLUDE_FROM_ALL

    algebra/curves/tests/test_groups.cpp
  )
  target_link_libraries(
    algebra_groups_test

    ff
    gtest_main
  )

  add_executable(
    algebra_fields_test
    EXCLUDE_FROM_ALL

    algebra/fields/tests/test_fields.cpp
  )
  target_link_libraries(
    algebra_fields_test

    ff
    gtest_main
  )

  add_executable(
    algebra_multiexp_test
    EXCLUDE_FROM_ALL

    algebra/scalar_multiplication/test_multiexp.cpp
  )
  target_link_libraries(
    algebra_multiexp_test

    ff
    gtest_main
  )

  include(CTest)
  add_test(
    NAME algebra_bilinearity_test
    COMMAND algebra_bilinearity_test
  )
  add_test(
    NAME algebra_groups_test
    COMMAND algebra_groups_test
  )
  add_test(
    NAME algebra_fields_test
    COMMAND algebra_fields_test
  )

  add_test(
    NAME algebra_multiexp_test
    COMMAND algebra_multiexp_test
  )

  add_dependencies(check algebra_bilinearity_test)
  add_dependencies(check algebra_groups_test)
  add_dependencies(check algebra_fields_test)
  add_dependencies(check algebra_multiexp_test)

  add_executable(
    multiexp_profile
    EXCLUDE_FROM_ALL

    algebra/scalar_multiplication/multiexp_profile.cpp
  )
  target_link_libraries(
    multiexp_profile

    ${OPENSSL_LIBRARIES}
    ff
  )

  add_executable(
    algebra_groups_profile
    EXCLUDE_FROM_ALL

    algebra/curves/profile/algebra_groups_profile.cpp
  )
  target_link_libraries(algebra_groups_profile ff)

  add_executable(
    algebra_groups_read_profile
    EXCLUDE_FROM_ALL

    algebra/curves/profile/algebra_groups_read_profile.cpp
  )
  target_link_libraries(algebra_groups_read_profile ff)

  add_dependencies(
    profile

    multiexp_profile
    algebra_groups_profile
    algebra_groups_read_profile
  )
endif()
