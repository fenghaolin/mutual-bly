# General Discription
1. This repository is a fork of [laseryuan/mutual-bly](https://github.com/laseryuan/mutual-bly), which is the code for the paper [Bensoussan, Yuan, and Liu (2023)](#BLY). The purpose of this fork is not to comment on its efficiency or to conduct further numerical experiments. Instead, it is to keep a copy of the version (commit 95c935b997576445fb6d3f9b12f2e1af4bd5494f) for easy future reference and to document the building and running environment of the code for reproducibility purposes. 
2. [Bensoussan, Yuan, and Liu (2023)](#BLY) is a paper published in *Operations Research*, a top journal in the field of operations research and management science. It studies the two-band impulse control of Brownian motion with the two-piece linear inventory cost structure ($f(x)=hx^++px^-$), which is a special case of the problem with a general convex $f(x)$ studied by [Feng and Muthuraman (2010)](#FM). [Bensoussan, Yuan, and Liu (2023)](#BLY) develop a new method called the splitting method and benchmark their method against the one proposed by [Feng and Muthuraman (2010)](#FM). The authors of [Bensoussan, Yuan, and Liu (2023)](#BLY) implement the splitting method and the benchmarking experiments in this codebase.
3. The source codes are the same as those in the [original repo](https://github.com/laseryuan/mutual-bly) (commit 95c935b997576445fb6d3f9b12f2e1af4bd5494f)
   - Only the README.md is modified to better document the steps to build and run the code.
   - A `.gitignore` file is added to ignore the built executable file `bin/Mutual-BLY`.
   - A `cpp_build.yml` file for Github Actions is added. Its path is `.github/workflows/cpp_build.yml`. This is added to attemp to run the code using Github Actions.

## Setup and Build: Memo by H. Feng
### Environment
- It has been tested on WSL 2 (Windows Subsystem Linux Version 2) with the Ubuntu 22.04 installed.
- It should also work in other Unix-like systems

### Dependencies
1. The compiler: clang
2. The GSL development libraries (libgsl-dev) 

- If any of them is missing, it is straightforwad to install the missing one. 
- For example, in a Debian-based system (e.g., Ubuntu), we can install libgsl-dev by:
```
sudo apt install libgsl-dev
``` 

- Note that Tthe Eigen library is also needed, but the files are already in 'include/eigen' of this repo.


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
clang++ -I ./include/eigen -I /usr/include/gsl -std=c++11 Main.cpp -o bin/Mutual-BLY /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgslcblas.so
```

## Note by H. Feng
- The code is written in C++ by the authors of [Bensoussan, Yuan, and Liu (2023)](#BLY). Without any modification of the source files, the project can be built successfully using the above `clang` command with warning message about with the style not following the requirements of ISO C++. However, if one wants to replace the clang compiler with the g++ one, it will fail with errors, some of which are easily fixable while others may not be. This note is a reminder for anyone interested in building and running the code.

## References
<a name="BLY"> [1] Bensoussan, A., Liu, J. J., & Yuan, J. (2023). A Splitting Method for Band Control of Brownian Motion: With Application to Mutual Reserve Optimization. Operations Research.</a>

<a name="FM"> [2] Feng, H., & Muthuraman, K. (2010). A computational method for stochastic impulse control problems. Mathematics of Operations Research, 35(4), 830-850.</a>

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
