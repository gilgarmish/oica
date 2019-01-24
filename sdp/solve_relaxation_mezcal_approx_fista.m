function D = solve_relaxation_mezcal_approx_fista(Fsbasis,G,mu,Dinit,maxiter,tolerance)
% Copyright: Anastasia Podosinnikova 2019

  [d,k] = size(Fsbasis);
  d = sqrt(d);

  D = Dinit;
  E = D;
  L = mu;
  t = 1;
  primal_vals = zeros(1,maxiter);
  dual_vals = zeros(1,maxiter);
  for iter = 1:maxiter
      temp = Fsbasis' * E(:);
      grad = -G + mu * reshape( Fsbasis * temp,d,d);
      E = E - (1/L) * ( grad );
      tnew = .5 * ( 1 + sqrt( 1 + 4 *t*t ) );
      [u,e] = eig((E+E')/2);
      u = real(u);
      e = real(diag(e));
      eproj = proj_simplex( e(:) );
      Dnew = u * diag( eproj ) * u'; D = (D+D')/2;
      E = Dnew + ( t - 1) / tnew * ( Dnew - D);
      D = Dnew;
      t = tnew;

      if mod(iter,10)==1
          temp = Fsbasis' * D(:);
          grad = -G + mu * reshape( Fsbasis * temp,d,d);
          primal_vals(iter) = -sum( G(:) .* D(:) ) + mu/2 * sum( temp.^2 );
          dual_vals(iter) = min(real(eig(grad))) - (mu/2) * sum( temp.^2);

          if ( (primal_vals(iter) - max(dual_vals(1:10:iter)) ) < tolerance )
              break;
          end
      end

  end

end

