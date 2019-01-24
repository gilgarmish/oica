function algname = get_algname( alg )
% Copyright: Anastasia Podosinnikova 2019
  algname = alg;
  
  if strcmp( alg, 'oica-clust' ), algname = 'OICA(QC)'; end
  if strcmp( alg, 'oica-ada' ), algname = 'OICA(QA)'; end
  if strcmp( alg, 'oica-semiada' ), algname = 'OverICA'; end
  if strcmp( alg, 'oica' ), algname = 'OverICA'; end
  
  if strcmp( alg, 'oica-quad-clust' ), algname = 'OICA(QC)'; end
  if strcmp( alg, 'oica-quad-clust-h' ), algname = 'OICA(QC-h)'; end
  if strcmp( alg, 'oica-quad-clust-km' ), algname = 'OICA(QC-km)'; end
  if strcmp( alg, 'oica-quad-ada' ), algname = 'OICA(QA)'; end
  if strcmp( alg, 'oica-quad-semiada' ), algname = 'OverICA(Q)'; end
  if strcmp( alg, 'oica-quad-semiada-h' ), algname = 'OverICA'; end %'OICA(QS)'; end
  if strcmp( alg, 'oica-quad-semiada-km' ), algname = 'OICA(QS-km)'; end
  if strcmp( alg, 'oica-gencov' ), algname = 'OICA(G)'; end
  if strcmp( alg, 'foobi' ), algname = 'FOOBI'; end
  if strcmp( alg, 'fpca' ), algname = 'Fourier PCA'; end
  if strcmp( alg, 'rica' ), algname = 'RICA'; end
  if strcmp( alg, 'rand' ), algname = 'RAND'; end
end
