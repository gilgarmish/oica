function D = majorize_minimize(G, Fs)
% Copyright: Anastasia Podosinnikova 2019
  p = sqrt( size(Fs,1) );
  Dinit = eye(p) / p;
  D = Dinit;
  mu = 5;
  maxiter = 100;
  tolerance = 1e-3;
  nmmmax = 100;
  
  iter = 1;
  while norm(D, 'fro') < 1
    u = extract_largest_eigenvector(G);
    Ginit = u*u';
    D = solve_relaxation_mezcal_approx_fista(Fs,Ginit,mu,Dinit,maxiter,tolerance);
    G = (D+D')/2;
    iter = iter + 1;
    if iter > nmmmax
      break;
    end
  end
end