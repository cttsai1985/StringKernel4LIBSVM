add_executable(stringkernel_tests
  test_kernels.cpp test_dataset.cpp test_options.cpp test_writer.cpp)
target_link_libraries(stringkernel_tests PRIVATE stringkernel_options Catch2::Catch2WithMain)
stringkernel_project_target(stringkernel_tests)
list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
include(Catch)
catch_discover_tests(stringkernel_tests DISCOVERY_MODE PRE_TEST)

add_test(NAME cli COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_cli.py"
  --converter "$<TARGET_FILE:StringKernelConverter>")
add_test(NAME libsvm_integration
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_libsvm_integration.py"
  --converter "$<TARGET_FILE:StringKernelConverter>"
  --trainer "$<TARGET_FILE:svm-train>" --predictor "$<TARGET_FILE:svm-predict>"
  --source "${PROJECT_SOURCE_DIR}" --cmake "${CMAKE_COMMAND}"
  --build "${PROJECT_BINARY_DIR}" --config "$<CONFIG>")
add_test(NAME separation_analysis
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_analysis.py"
  --converter "$<TARGET_FILE:StringKernelConverter>"
  --trainer "$<TARGET_FILE:svm-train>" --predictor "$<TARGET_FILE:svm-predict>"
  --source "${PROJECT_SOURCE_DIR}")
add_test(NAME dependency_pins
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_dependency_pins.py")
add_test(NAME dependency_bootstrap
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_fetch_libsvm.py")
set_tests_properties(cli libsvm_integration separation_analysis dependency_pins dependency_bootstrap
  PROPERTIES TIMEOUT 120)
