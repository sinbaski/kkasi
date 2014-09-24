clear all
close all
q = 0.1;

for sig = [0.1, 0.5, 1]
    load(sprintf(['../matfys/data/sv/normal_ret/lognormal_vol/Eig-' ...
                  'sig%.4f-phi%.4f.mat'],sig, phi));

    ev = reshape(ev, 1, prod(size(ev)));
    
end



