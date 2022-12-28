function [interps,coefs] = getTargetFunctionExt(rho1,rho2,kernelType,targetHyper)
% Generating the target function through a kernel matrix based on the norm of samples in 3 dimensional space 
%% Generate interpolating inputs for target function

n1 = 16;
n2 = 8;
theta1Range = linspace(0,2*pi,n1);
theta2Range = linspace(0,2*pi,n2);
m = n1*n2;
[theta1Interps,theta2Interps] = meshgrid(theta1Range,theta2Range);


u = (rho1+rho2*cos(theta2Interps)).*cos(theta1Interps);
v = (rho1+rho2*cos(theta2Interps)).*sin(theta1Interps);
w = rho2*sin(theta2Interps);
%% Generate interpolating values;
p = @(x,y,z) 1/8*(x^5-10*x^3*y^2+5*x*y^4)*(x^2+y^2-60*z^2);
p_samps = zeros(m,1);
figure()
plot3(u,v,w,'ko')

%% making the list of centers
interps = zeros(m,3);
index = 1;
for ii = 1:n2
    for jj = 1:n1
        interps(index,:) = [u(ii,jj),v(ii,jj),w(ii,jj)];
        p_samps(index) = p(u(ii,jj),v(ii,jj),w(ii,jj));
        index = index + 1;
    end
end

%% calculating coefficients for target function
A = zeros(m,m);
for pp = 1:m
    for nn = 1:m
        A(pp,nn) = kernel(kernelType,interps(nn,:),interps(pp,:),targetHyper);
    end
end
title(cond(A))
coefs = pinv(A)*p_samps;
