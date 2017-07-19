StringKernelConverter.exe -train data.txt -k 5 > libsvm_train.txt
StringKernelConverter.exe -train data.txt -test test_data.txt -k 5 > libsvm_test.txt