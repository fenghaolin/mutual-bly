# This is forked from 
- The source codes are the same as the [original repo](https://github.com/laseryuan/mutual-bly) (commit 95c935b997576445fb6d3f9b12f2e1af4bd5494f)
- Only the README.md is modified to better document the steps to build and run it.
- A `.gitignore` file is added to ignore the built executable file `bin/Mutual-BLY`
- A `cpp_build.yml` file for Github Actions is added. Its path is `.github/workflows/cpp_build.yml`
	- This is added to attemp to run the code using Github Actions

## Setup and Build: Memo by H. Feng
### Environment
- It has been tested on WSL 2 (Windows Subsystem Linux Version 2) with the Ubuntu 22.04 installed.
- It should also work in other Unix-like systems

### Dependencies
1. The compiler clang
2. The GSL development libraries (libgsl-dev) 

- If any of them is missing, it is straightforwad to install the missing one. 
- For example, in a Debian-based system (e.g., Ubuntu), we can install libgsl-dev by:
```
sudo apt install libgsl-dev
``` 

- Note that Tthe Eigen library is also needed, but the files are already in 'include/eigen' of this repo


### Building and running it
- Assuming 
	1. The headers of the 'gsl' library is in the directory '/usr/include/gsl'
	2. The shared libraries 'libgsl.so' and 'libgslcblas.so' for the GSL library can be found in the directory '/usr/lib/x86_64-linux-gnu/'
		- If this is not the case, one can find them using the following two commands, repectively,
		```
		sudo find /usr/lib/ /usr/local/lib/ /lib/ -name libgsl.so 2>/dev/null
		sudo find /usr/lib/ /usr/local/lib/ /lib/ -name libgslcblas.so 2>/dev/null
		```
- Step 1. git clone the repo and change directory to the local folder using
```
git clone https://github.com/fenghaolin/mutual-bly.git mutual-bly-OR2023 && cd mutual-bly-OR2023
```

- Step 2. Buildi the executable file `Mutual-BLY` and put it into `bin`:
```
clang++ -I ./include/eigen -I /usr/include/gsl \
    -std=c++11 Main.cpp -o bin/Mutual-BLY \ 
    /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgslcblas.so
```

# In below is the README.md of the [original](https://github.com/laseryuan/mutual-bly) repo 
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
