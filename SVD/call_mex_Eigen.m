% Make sure to modify and point to the correct location of your Eigen
% installation.
mex CXXFLAGS="\$CXXFLAGS -std=c++11 -fopenmp" CXXOPTIMFLAGS='\$CXXOPTIMFLAGS -Ofast -DNDEBUG' LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp -O2" -lgomp -I"./Eigen"  batchSVD3x3Eigen.cpp
