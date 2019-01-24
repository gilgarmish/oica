function Ds = Ds_from_ds(ds)
% Copyright: Anastasia Podosinnikova 2019
  [p,k] = size(ds);
  Ds = zeros(p^2, k);
  for i=1:k
    Ds(:,i) = vec( ds(:,i)*ds(:,i)' );
  end
end