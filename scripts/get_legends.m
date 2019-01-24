function legends = get_legends( algs )
% Copyright: Anastasia Podosinnikova 2019
  nalg = length(algs);
  legends = cell(nalg,1);
  for i = 1:nalg
    legends{i} = get_algname( algs{i} );
  end
end