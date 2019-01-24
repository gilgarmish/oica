function [u,e] = extract_largest_eigenvector(D)
% Copyright: Anastasia Podosinnikova 2019
  D = (D+D')/2;
  [u,e] = eig(D);
  %[u,e] = eigs(D,1);
  [a,b] = max(abs(diag(e)));
  u = real(u(:,b(1)));
  e = a(1);
end
