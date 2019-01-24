function C = gencov(X, omega)
% Copyright: Anastasia Podosinnikova 2019

  if numel(omega)==1
    p = size(X,1);
    omega = omega*ones(p,1)/p;
  end

  n = size(X,2);
  
  proj = X'*omega;
  eproj = exp(proj);
  % genexp
  Eomega = (X*eproj) / sum(eproj); 
  
  C = X*sparse(1:n,1:n,eproj)*X' ;
  C = C / sum(eproj);
  C = C - Eomega * Eomega';
  C = C(:);
end
