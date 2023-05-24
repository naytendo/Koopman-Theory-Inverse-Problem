function [interps,coefs] = getTargetFunctionUniform(rho1,rho2,kernelType,targetHyper,nu)
%% Generate interpolating inputs for target function
m = 100;
interps = zeros(m,3);
p_samps = zeros(m,1);
ii = 1;

p = @(x,y,z) 1/8*(x^5-10*x^3*y^2+5*x*y^4)*(x^2+y^2-60*z^2);

gmax = rho2 * (rho2 + rho1);

while ii <= m
    theta1 = 2*pi*rand;
    theta2 = 2*pi*rand;
    x = (rho1+rho2*cos(theta2))*cos(theta1);
    y = (rho1+rho2*cos(theta2))*sin(theta1);
    z = rho2*sin(theta2);
    g = rho2*abs(rho1+cos(theta2)*rho1);

    if g >= gmax*rand
        interps(ii,:) = [x,y,z];
        ii = ii +1;
        p_samps(ii) = p(x,y,z);
    end
    
end



%% calculating coefficients for target function
A = zeros(m,m);
for pp = 1:m
    for nn = 1:m
        A(pp,nn) = kernel(kernelType,interps(nn,:),interps(pp,:),targetHyper,nu);
    end
end
% title(cond(A))
coefs = pinv(A)*p_samps;

