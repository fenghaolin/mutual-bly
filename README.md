## Setup
```
git clone https://github.com/laseryuan/mutual-bly && cd mutual-bly  
docker run -it -v $PWD:/home/app/app lasery/cling bash
```

## Run
It took about 3 minutes to finish all the test cases on Intel Core i7-2600 @ 3.40GHz.

The performance comparison between our method and the benchmark is recorded in
the benchmark.txt file.
```
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib bin/Mutual-BLY
cat output/benchmark.txt
```

## Build
```
clang++ -I /home/app/app/include/eigen -I /usr/local/include/gsl -std=c++11 Main.cpp -o bin/Mutual-BLY /usr/local/lib/libgsl.so /usr/local/lib/libgslcblas.so
```
