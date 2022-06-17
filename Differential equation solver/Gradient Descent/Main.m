A_hom = cell2mat(struct2cell(load('A_hom.mat')));
sol_hom = Solver(A_hom);

A_sm = cell2mat(struct2cell(load('A_sm.mat')));
sol_sm = Solver(A_sm);

A_sf = cell2mat(struct2cell(load('A_sf.mat')));
sol_sf = Solver(A_sf);