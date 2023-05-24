
function [interps,coefs] = getTargetFunctionNew(rho1,rho2,targetType,targetScale,targetPar)
%% Generate interpolating inputs for target function
n1 = 17;
n2 = 9;

theta1Range = linspace(0,2*pi,n1);
theta2Range = linspace(0,2*pi,n2);
theta1Range = theta1Range(1:end-1);
theta2Range = theta2Range(1:end-1);
m = (n1-1)*(n2-1);
[theta1Interps,theta2Interps] = meshgrid(theta1Range,theta2Range);

u = (rho1+rho2*cos(theta2Interps)).*cos(theta1Interps);
v = (rho1+rho2*cos(theta2Interps)).*sin(theta1Interps);
w = rho2*sin(theta2Interps);
%% Generate interpolating values;
p = @(x,y,z) 1/8*(x^5-10*x^3*y^2+5*x*y^4)*(x^2+y^2-60*z^2);
p_samps = zeros(m,1);
% figure()
% plot3(u,v,w,'ko','markerFaceColor','k')
% hold on

%% making the list of centers
interps = zeros(m,3);
index = 1;
for ii = 1:n2-1
    for jj = 1:n1-1
        interps(index,:) = [u(ii,jj),v(ii,jj),w(ii,jj)];
        p_samps(index) = p(u(ii,jj),v(ii,jj),w(ii,jj));
        index = index + 1;
    end
end

%% calculating coefficients for target function
A = real(kermat(interps,interps,targetType,targetPar,targetScale));
% title(cond(A))
coefs = A\p_samps;

