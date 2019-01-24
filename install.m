function install( ismex )

% add dependencies

addpath(genpath(pwd))


% compile mex files

if nargin==1 && ismex==1

  cd comparison
    cd minfunc
      mexAll
    cd ..
    cd ojd
      mex ojd_in.cpp -largeArrayDims;
    cd ..
  cd ..

  cd cumulants
    mex quadricov_in.cpp -largeArrayDims;
  cd ..

end

end
