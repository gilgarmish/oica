function [ds, eigmaxes] = approx_ds_from_Ds(Ds)
% Copyright: Anastasia Podosinnikova 2019
  k = size(Ds,2);
  p = sqrt(size(Ds,1));
  ds = zeros(p,k);
  eigmaxes = zeros(k,1);
  
  for i = 1:k
    D = reshape( Ds(:,i),p,p ); 
    [u,e] = extract_largest_eigenvector(D);
    eigmaxes(i) = e;
    ds(:,i) = u;
  end
  
end