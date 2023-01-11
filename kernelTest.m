targetType = 'matern';
targetHyper = 1;

center = [0;0];
x = [zeros(1,1000),linspace(0,5,1000)];
b = 9;
dim = 2;
nu = 0.5*(b-dim/2);

kern = zeros(1,length(1));
for ii = 1:length(x)
    kern(ii) = kernel(targetType,center,x(:,ii),targetHyper,nu);
end

plot(x,kern)