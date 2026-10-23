# Sources are prepared explicitly; production configuration never downloads or
# invokes Python. Validate the exact bytes on every configuration, including
# offline builds and locally modified dependency caches.
set(libsvm_dir "${PROJECT_SOURCE_DIR}/third_party/libsvm")
file(READ "${PROJECT_SOURCE_DIR}/cmake/libsvm.json" libsvm_manifest)
string(JSON libsvm_file_count LENGTH "${libsvm_manifest}" files)
math(EXPR libsvm_last_file "${libsvm_file_count} - 1")
foreach(index RANGE ${libsvm_last_file})
  string(JSON name MEMBER "${libsvm_manifest}" files ${index})
  string(JSON expected_hash GET "${libsvm_manifest}" files "${name}" sha256)
  if(NOT EXISTS "${libsvm_dir}/${name}")
    message(FATAL_ERROR
      "Missing LIBSVM dependency: ${name}. Run python3 scripts/fetch_libsvm.py "
      "from the repository root (Windows: python scripts/fetch_libsvm.py).")
  endif()
  file(SHA256 "${libsvm_dir}/${name}" actual_hash)
  if(NOT actual_hash STREQUAL expected_hash)
    message(FATAL_ERROR
      "LIBSVM SHA-256 mismatch: ${name}. Remove the modified file and run "
      "python3 scripts/fetch_libsvm.py (Windows: python scripts/fetch_libsvm.py).")
  endif()
endforeach()
