function [Qinit, kmax, eps, isdebug] = initialize_defaults(m)
% Copyright: Anastasia Podosinnikova 2019
  Qinit = eye(m);
  kmax = 1000;
  eps = 1e-5;
  isdebug = 0;
end