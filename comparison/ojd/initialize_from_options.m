function [Qinit, kmax, eps, isdebug] = initialize_from_options( options, m )
% Copyright: Anastasia Podosinnikova 2019

  if ~isstruct( options ), error('wrong input'), end

  if isfield( options, 'Qinit' )
    Qinit = options.Qinit;
  else
    Qinit = eye(m);
  end

  if isfield( options, 'kmax' )
    kmax = options.kmax;
  else
    kmax = 1000;
  end

  if isfield( options, 'eps' )
    eps = options.eps;
  else
    eps = 1e-8;
  end

  if isfield( options, 'isdebug' )
    isdebug = options.isdebug;
  else
    isdebug = 0;
  end

end
