function reproduce_fig_1_phase_transition
% This script reproduces the phase transition plots (Figure 1)
% Dependency: cvx optimization toolbox (http://cvxr.com/cvx/)
% This script as is, without parallelization, will run for many hours

% Copyright: Anastasia Podosinnikova 2019

  ps = 10:5:50;
  ks = 10:10:950;
  ptvals = zeros( length(ps), length(ks) );
  nruns = 10;
  
  ks4p = cell( length(ps), 1 );
  ks4p{1} = 10:10:50;
  ks4p{2} = 10:10:120;
  ks4p{3} = 10:10:210;
  ks4p{4} = 10:10:220;
  ks4p{5} = 10:10:400;
  ks4p{6} = 10:10:500;
  ks4p{7} = 10:10:650;
  ks4p{8} = 10:10:800;
  ks4p{9} = 10:10:950;
  
  for i = 1:length(ps)
    p = ps(i);
    ksloc = ks4p{i};
    
    filepath = strcat( pwd, '/expres/pt/p', num2str(p),'.mat' );
    if exist( filepath, 'file' ) == 2
      rr = load(filepath);
      successes = rr.successes;
    else
      successes = compute(p, ksloc, nruns);
      save( filepath, 'successes' )
    end
    
    for j = 1:length(ksloc)
      ptvals( logical(ps==p), logical(ks==ksloc(j)) ) = successes(j);
    end
  end
  
  make_plot(ptvals)

end



function successes = compute(p, ks, nruns)
  n = length(ks);
  successes = zeros(n,1);
  
  for i = 1:n
    k = ks(i);
    
    ds = sample_mixing_matrix(p,k);
    Ds = zeros(p^2,k);
    for j = 1:k
      Ds(:,j) = vec(ds(:,j)*ds(:,j)');
    end
    
    for irun = 1:nruns

      u = randn(p,1); u = u/norm(u);
      G = u*u';
      traces = Ds' * vec(G);

      cvx_begin quiet
      expression B(p,p) 
      variable q(k,1)
      maximize ( traces' * q  )
      subject to
        B = zeros(p,p);
        for j=1:k
          B = B + q(j) * reshape( Ds(:,j), p,p );
        end
        trace(B) == 1
        B == semidefinite(p)
      cvx_end
      
      if abs( max(Ds'*vec(B)) - 1 ) < 0.001
        successes(i) = successes(i) + 1;
      end
      
    end
    
    successes(i) = successes(i) / nruns;
  end
end


function make_plot(ptvals)

  ff=figure; hold on
    
    screensize = get( groot, 'Screensize' );
    position = [0 0 screensize(3) screensize(4)];
    set(ff, 'Position', position)
    
    % font type for text and ticks
    fontname = 'Times New Roman';
    set(gca, 'Fontname', fontname)
    
    % font size for title, axis, ticks
    fontsize = 60;
    set(gca, 'FontSize', fontsize)
    
    pcolor( ptvals ), colorbar, colormap(flipud(gray))
    shading flat
    
    xlim([1 95])
    ylim([1 9])
    
    yys = .5:1:9.5;

    ys = [0 10:5:50];

    xxs = (ys .* (ys+1) / 2) / 10;

    zzs = (ys.^2 / 4) / 10;

    plot( xxs, yys, 'Color', 'b', 'LineWidth', 5 )
    plot( zzs, yys, 'Color', 'r', 'LineWidth', 5 )
    
    set(gca,'YTick', 1:9)
    set(gca,'YTickLabel', 10:5:50)
    
    set(gca,'XTick', 10:10:90)
    set(gca,'XTickLabel', 100:100:900)

    xlabel('Latent Dimension ($k$)', 'Interpreter', 'latex')
    ylabel('Observed Dimension ($p$)', 'Interpreter', 'latex')
    
    box on

end
