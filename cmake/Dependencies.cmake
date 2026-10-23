include(FetchContent)

# Test-only dependency. For offline tests set FETCHCONTENT_SOURCE_DIR_CATCH2
# to a previously downloaded, verified Catch2 v3.16.1 source directory.
set(CATCH_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
set(CATCH_INSTALL_EXTRAS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(catch2
  URL https://codeload.github.com/catchorg/Catch2/tar.gz/refs/tags/v3.16.1
  URL_HASH SHA256=d16bd1b8b2364918bd472cb216394acb93bde6784890df0d91bccfb3f28d777d
  DOWNLOAD_EXTRACT_TIMESTAMP FALSE)
FetchContent_MakeAvailable(catch2)
# Catch2's compiled feature support must match the C++17 test headers on every
# compiler, including AppleClang whose default standard can be older.
target_compile_features(Catch2 PUBLIC cxx_std_17)
