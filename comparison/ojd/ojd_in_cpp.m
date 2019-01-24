function [Q, diags] = ojd_in_cpp(A, options)
% a wrapper to call ojod_in.ccp
% ojod : orthogonal joint off diagonalization
% see nojod_in_cpp.m for the description of options

% Copyright: Anastasia Podosinnikova 2019

  % process input
  if (nargin < 1) || (nargin > 2)
    error('Wrong input')
  end
  [m, nm] = size(A); n = nm/m;
  if nargin==1
    [Qinit, kmax, eps, isdebug] = initialize_defaults(m);
  else
    [Qinit, kmax, eps, isdebug] = initialize_from_options( options, m );
  end
  
  % TODO: initializing with a diagonalizer of one of the matrices
  
  % main
  B = A;
  for i=1:n
    inds = (i-1)*m+1 : i*m;
    B(:,inds) = Qinit'*B(:,inds)*Qinit;
  end
  
  [out1, out2] = ojd_in( B(:), m, n, eps, Qinit(:), kmax, isdebug); 
  
  [Q, diags] = reshape_the_output_of_jd_in_cpp(out1, out2, m, n);
  
end

function [Q, diags] = reshape_the_output_of_jd_in_cpp(out1, out2, m, n)
% reshape the output of jd_in.cpp
  Q = reshape(out1, m, m);
  diags_temp = reshape(out2, m, m*n);
  diags = cell(1,n);
  for i = 1:n, 
    diags{i} = diags_temp(:, (i-1)*m+1 : i*m); 
  end
end
