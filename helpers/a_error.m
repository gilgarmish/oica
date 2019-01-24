function [perf, perfout] = a_error(ds_est, ds)
% Copyright: Anastasia Podosinnikova 2019
  perf_type = 4;
  [perf, perfout] = evaluation_perf(ds_est, ds, perf_type);
end