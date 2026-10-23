@echo off
setlocal
pushd "%~dp0" || exit /b 1

set "BIN=%~1"
if not defined BIN set "BIN=build\release\Release"
if not exist "%BIN%\StringKernelConverter.exe" if "%~1"=="" set "BIN=build\release"
if not exist "%BIN%\StringKernelConverter.exe" (
  echo Converter not found. Pass the directory containing the built executables. 1>&2
  goto :fail
)
if not exist "build\example" mkdir "build\example"
if not exist "build\example" goto :fail

"%BIN%\StringKernelConverter.exe" -train examples\regression\train_sequences.tsv -k 5 > "build\example\train.kernel"
if errorlevel 1 goto :fail
"%BIN%\StringKernelConverter.exe" -train examples\regression\train_sequences.tsv -test examples\regression\test_sequences.tsv -k 5 > "build\example\test.kernel"
if errorlevel 1 goto :fail
"%BIN%\svm-train.exe" -s 3 -t 4 "build\example\train.kernel" "build\example\model"
if errorlevel 1 goto :fail
"%BIN%\svm-predict.exe" "build\example\test.kernel" "build\example\model" "build\example\predictions.txt"
if errorlevel 1 goto :fail

popd
exit /b 0

:fail
popd
exit /b 1
