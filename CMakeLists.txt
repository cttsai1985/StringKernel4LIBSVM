cmake_minimum_required(VERSION 3.25)
project(StringKernel4LIBSVM VERSION 1.0.0 LANGUAGES CXX)

include(CTest)
include(GNUInstallDirs)

option(STRINGKERNEL_SANITIZERS "Enable address and undefined-behavior sanitizers" OFF)
option(STRINGKERNEL_COVERAGE "Enable GCC coverage for project sources" OFF)

if(STRINGKERNEL_SANITIZERS AND (WIN32 OR NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang"))
  message(FATAL_ERROR "Sanitizers require GCC or Clang on a supported Unix platform")
endif()
if(STRINGKERNEL_COVERAGE AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR "Coverage requires GCC")
endif()
if(STRINGKERNEL_SANITIZERS AND STRINGKERNEL_COVERAGE)
  message(FATAL_ERROR "Use separate builds for sanitizers and coverage")
endif()

function(stringkernel_project_target target)
  target_compile_features(${target} PRIVATE cxx_std_17)
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX /permissive- /utf-8)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
  if(STRINGKERNEL_SANITIZERS)
    target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(${target} PUBLIC -fsanitize=address,undefined)
  endif()
  if(STRINGKERNEL_COVERAGE)
    target_compile_options(${target} PRIVATE --coverage -O0 -g)
    target_link_options(${target} PUBLIC --coverage)
  endif()
endfunction()

add_library(stringkernel STATIC StringKernels.cpp Dataset.cpp LibsvmWriter.cpp)
target_compile_features(stringkernel PUBLIC cxx_std_17)
target_include_directories(stringkernel PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/StringKernel4LIBSVM>)
stringkernel_project_target(stringkernel)

add_library(stringkernel_options STATIC CliOptions.cpp)
target_link_libraries(stringkernel_options PUBLIC stringkernel)
stringkernel_project_target(stringkernel_options)

include(cmake/Libsvm.cmake)
add_library(libsvm STATIC "${libsvm_dir}/svm.cpp")
target_include_directories(libsvm PUBLIC "${libsvm_dir}")
set_target_properties(libsvm PROPERTIES CXX_EXTENSIONS OFF)
set_source_files_properties("${libsvm_dir}/svm-train.c" "${libsvm_dir}/svm-predict.c"
  PROPERTIES LANGUAGE CXX)
add_executable(svm-train "${libsvm_dir}/svm-train.c")
add_executable(svm-predict "${libsvm_dir}/svm-predict.c")
target_link_libraries(svm-train PRIVATE libsvm)
target_link_libraries(svm-predict PRIVATE libsvm)
if(MSVC)
  target_compile_definitions(libsvm PRIVATE _CRT_SECURE_NO_WARNINGS)
  target_compile_definitions(svm-train PRIVATE _CRT_SECURE_NO_WARNINGS)
  target_compile_definitions(svm-predict PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

add_executable(StringKernelConverter StringKernelConverter.cpp)
target_link_libraries(StringKernelConverter PRIVATE stringkernel_options)
target_include_directories(StringKernelConverter PRIVATE "${libsvm_dir}")
target_compile_definitions(StringKernelConverter PRIVATE STRINGKERNEL_VERSION="${PROJECT_VERSION}")
stringkernel_project_target(StringKernelConverter)

install(TARGETS stringkernel EXPORT stringkernelTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS StringKernelConverter svm-train svm-predict
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES StringKernels.h Dataset.h LibsvmWriter.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/StringKernel4LIBSVM)
install(EXPORT stringkernelTargets NAMESPACE stringkernel::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/StringKernel4LIBSVM)
install(FILES THIRD_PARTY_NOTICES.md "${libsvm_dir}/COPYRIGHT"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/StringKernel4LIBSVM)

if(BUILD_TESTING)
  include(cmake/Dependencies.cmake)
  find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
  add_subdirectory(tests)
endif()
