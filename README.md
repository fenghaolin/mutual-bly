## Development enviorment
git clone https://github.com/laseryuan/mutual-bly && cd mutual-bly
docker run -it -v $PWD:/home/app/app lasery/cling bash

## Build
clang++ -I /home/app/app/include/eigen -I /usr/local/include/gsl -std=c++11 Main.cpp -o bin/Mutual-BLY /usr/local/lib/libgsl.so /usr/local/lib/libgslcblas.so

## Run
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib bin/Mutual-BLY
