function reproduce_fig_2_3_6_7_synthetic_exps
% This script reproduces experiments with synthetic data (Figures 2,3,6,7) 
% Dependency: FOOBI algorithm; to obtain the code, please contact 
%             Lieven De Lathauwer (http://homes.esat.kuleuven.be/~delathau/)
% This script as is, whithout parallelization, will run for a few hours

% Copyright: Anastasia Podosinnikova 2019

  % FIGURE 2
  randn('state',0); %#ok
  rand('state',0);  %#ok
  p = 10;
  ks = 5:5:60;
  nrep = 10;
  algs = {'rand', 'foobi', 'oica-semiada'};
  experiment_asymptotic( p, ks, nrep, algs )


  % FIGURE 3 (first two plots)
  randn('state',0); %#ok
  rand('state',0);  %#ok
  p = 15;
  k = 30;
  nrep = 10;
  algs = { 'foobi', 'oica', 'fpca', 'rand' };
  ns = 1e3 : 1e3 : 1e4;
  isshort = 1;
  experiment_fixedk( p, k, nrep, algs, ns, isshort )

  % FIGURES 3 AND 7
  randn('state',0); %#ok
  rand('state',0);  %#ok
  ns = 1e4 : 1e4 : 2e5;
  algs = { 'foobi', 'oica', 'fpca', 'oica-quad-semiada' };
  isshort = 0;
  experiment_fixedk( p, k, nrep, algs, ns, isshort )


  % FIGURE 4 - RUNTIME
  randn('state',0); %#ok
  rand('state',0);  %#ok
  p = 20;
  ks = 20:20:100;
  n = 1e5;
  nrep = 10;
  algs = {'foobi', 'fpca', 'oica'};
  experiment_fixedn( p, ks, nrep, algs, n );

end