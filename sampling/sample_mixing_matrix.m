function ds = sample_mixing_matrix(p,k)
% Copyright: Anastasia Podosinnikova 2019
  ds = randn(p,k);
  for i=1:k, ds(:,i) = ds(:,i)/norm(ds(:,i),2); end
end