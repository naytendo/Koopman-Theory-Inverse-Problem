
function [interps,coefs] = getTargetFunctionAlongOrbit2(rho1,rho2,kernelType,targetScale,RBFpar)
%% Generate interpolating inputs for target function
n1 = 1000;

theta = linspace(0,2*pi,n1);


u = (rho1+rho2*cos(6*theta)).*cos(theta);
v = (rho1+rho2*cos(6*theta)).*sin(theta);
w = rho2*sin(6*theta);
xOrbit = [u',v',w']; 
%% Generate interpolating values;
p = @(x,y,z) 1/8*(x^5-10*x^3*y^2+5*x*y^4)*(x^2+y^2-60*z^2);
% figure()
% plot3(u,v,w,'ko','markerFaceColor','k')
% hold on
p_Plot = zeros(length(u),1);
for ss = 1:length(u)
    p_Plot(ss) = p(u(ss),v(ss),w(ss)); 
end
%% making the list of interps
interps = zeros(length(xOrbit),3);
p_samps = zeros(length(xOrbit),1);
Np = 1;
interps(Np,:) = xOrbit(1,:);
separation = 0.4;
for ii = 1:length(xOrbit)
    check = 0;
    for jj = 1:Np
        if norm(interps(jj,:)-xOrbit(ii,:)) > 2*separation
            check = check +1;
        end
    end
    if check == Np
        interps(Np+1,:) = xOrbit(ii,:);
        p_samps(Np) = p(u(ii),v(ii),w(ii));
        Np = Np+1;
    end
end

p_samps = p_samps(1:Np);
interps = interps(1:Np,:);

%% calculating coefficients for target function
A = kermat(interps,interps,kernelType,RBFpar,targetScale);
% title(cond(A))
coefs = A\p_samps;

