function Ds_est = cluster_Dss(Dss, k, ctype)
% Copyright: Anastasia Podosinnikova 2019

  % since the atoms are scaling (hence, sign) invariant
  % we first align them all to point to the same direction
  
  D1 = Dss(:,1);
  for i = 2:size(Dss,2)
    Di = Dss(:,i);
    Dss(:,i) = sign( D1'*Di) * Di;
  end
  
  Ds_est = extract_clusters(Dss, k, ctype);
  
end

function Ds_temp = extract_clusters(DD, nclust, ctype)
  p = sqrt(size(DD,1));
  
  if strcmp(ctype,'h') % hierarchical clustering
    cc = clusterdata(DD', nclust);
  end
  if strcmp(ctype,'km') % k-means++
    cc = kmeans(DD, nclust);
  end
  
  Ds_temp = zeros(p^2, nclust);
  for i = 1:nclust
    DDi = DD(:, logical(cc==i));
    Ds_temp(:,i) = DDi(:,1); %mean(DDi,2);
  end
end
