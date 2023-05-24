function mat=kermat(X, Y,RBFtype,RBFPar,RBFscale)
% creates kernel matrix
% for two point sets X and Y

[~, dx]=size(X);
[~, dy]=size(Y);
if dx~=dy
    error('Unequal space dimension for kermat arguments');
end
mat=frbf(distsqh(X,Y)./RBFscale^2,0,RBFtype,RBFPar);
end