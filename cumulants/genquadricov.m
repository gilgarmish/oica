function Q = genquadricov(X, u)
%************************************************%%%%%%%%%%%%%%%%%%%%%%%%%%
% X: p x N, p dimension, N number of samples
% u: p, processing point
% Q: p*p x p*p, matricized 4-th order generalized cumulant
%
% The matricization rule:
% Q( (i1-1)*p + i2, (i3-1)*p + i4 ) = CUM(i1,i2,i3,i4);
% (i.e. row-wise for the the first two and the last two dimensions)
%
% THE OUTPUT IS A COMPLEX MATRIX
%************************************************%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright: Anastasia Podosinnikova 2019

  [p,N] = size(X);
  
  Xu = X'*u;
  expXu = exp(1i*Xu);
  M0 = sum( expXu ) / N;
  
  M1 = (X * expXu ) / (N*M0);
  
  XC = X - repmat( M1, 1, N ); % "generalized zero-mean"
  
  M2 = zeros(p);   % "generalized order-2 moment"
  M4 = zeros(p*p); % "generalized order-4 moment"
  temp = zeros(p);
  
  for n = 1:N
    xn = XC(:,n);
    temp = xn*conj(ctranspose(xn)); % this is correct; see comment below
    M2 = M2 + expXu(n)*temp;
    M4 = M4 + (expXu(n)*temp(:))*conj(ctranspose(temp(:))); 
    % !!!!! The transpose of a complex number is the conjugate transpose in
    % matlab. We want just transpose !!!!!
  end
  
  M2 = M2 / (N*M0); M4 = M4 / (N*M0);
  Q = zeros(p*p);
  
  
  for i3 = 1:p
    for i4 = 1:p
      icol = (i3-1)*p+i4;
      temp(:) = M4(:,icol);
      % temp = temp - M2(i3,i4)*M2' - M2(:,i3)*M2(:,i4)' - M2(:,i4)*M2(:,i3)';
      temp = temp - M2(i3,i4)*M2 - M2(:,i3)*M2(i4,:) - M2(:,i4)*M2(i3,:); % since M2 must be symmetric!
      Q(:,icol) = temp(:);
    end
  end
  
end
