function [Esbasis, Fsbasis] = extract_basis(Es, k)
% Copyright: Anastasia Podosinnikova 2019
  [Q,~] = qr(Es);
  Esbasis = Q(:,1:k);
  Fsbasis = Q(:,k+1:end);
end