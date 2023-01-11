rho1 = 1;
rho2 = 1/3;
targetType = 'matern52';
targetHyper = 0.5;
[interps,coefs] = getTargetFunctionExt(rho1,rho2,targetType,targetHyper);

%% Generating target function surface
% Note: target function at a point is calculated by coefficients times the kernel vector 
meshRes = 0.1;
theta1Range = 0:meshRes:2*pi+meshRes;
theta2Range = 0:meshRes:2*pi+meshRes;

[theta1Mesh,theta2Mesh] = meshgrid(theta1Range,theta2Range);

u_Mesh = (rho1+rho2*cos(theta2Mesh)).*cos(theta1Mesh);
v_Mesh = (rho1+rho2*cos(theta2Mesh)).*sin(theta1Mesh);
w_Mesh = -rho2*sin(theta2Mesh);



targetFunc_Mesh = zeros(size(u_Mesh));
% 
for ii = 1:size(targetFunc_Mesh,1)
    for jj = 1:size(targetFunc_Mesh,2)
    kernVector = zeros(length(interps),1);
        for cc = 1:length(interps)
            kernVector(cc) = kernel(targetType,interps(cc,:),[u_Mesh(ii,jj),v_Mesh(ii,jj),w_Mesh(ii,jj)],targetHyper);
        end
    targetFunc_Mesh(ii,jj) = coefs'*kernVector;
    end
end
surf(u_Mesh,v_Mesh,w_Mesh,targetFunc_Mesh,'EdgeColor','none','FaceAlpha',0.9);
colorbar
title('target Function with interpolating points')


figure()

%%
n = 2;

gamma = 0.85;
delta = 0.05;
q = sqrt((0.5/delta)^2-1);
alpha0 = [q
1];
tau = 1.25;

[boolean,kProb,checks] = isDiophantine(delta,alpha0,gamma,tau);

alpha = alpha0/norm(alpha0);
if boolean
%     T = (1+n^2*factorial(n))^(tau+1)/gamma/delta;
    T = sqrt(q^2+1);
    t = 1:0.01:T;
    theta1_Orb = alpha(1)*t;
    theta2_Orb = alpha(2)*t;
    mod_theta1_Orb = zeros(length(theta1_Orb),1);
    mod_theta2_Orb = zeros(length(theta2_Orb),1);
    for ii = 1:length(theta1_Orb)
        mod_theta1_Orb(ii) = mod(theta1_Orb(ii),1);
        mod_theta2_Orb(ii) = mod(theta2_Orb(ii),1);
    end
    

    w1 = alpha(1);
    w2 = alpha(2);
    x_Orb = (rho1+rho2*cos(2*pi*w2*t)).*cos(2*pi*w1*t);
    y_Orb = (rho1+rho2*cos(2*pi*w2*t)).*sin(2*pi*w1*t);
    z_Orb = -rho2*sin(2*pi*w2*t);
    plot3(x_Orb,y_Orb,z_Orb,'k');
    hold on
    surf(u_Mesh,v_Mesh,w_Mesh,targetFunc_Mesh,'EdgeColor','none','FaceAlpha',0.9);
    colorbar
    title('\delta filling orbit on Torus in 3D')

    figure()
    plot(mod_theta1_Orb,mod_theta2_Orb,'.','markersize',0.5)
    grid on
    xlabel('\theta_1')
    ylabel('\theta_1')
    title('\delta filling orbit in 2D intrinsic coordinates')
end



figure()
surf(theta1Mesh/2/pi,theta2Mesh/2/pi,targetFunc_Mesh,'EdgeColor','none')
title('Target Function over 2D Torus')
xlabel('\theta_1')
ylabel('\theta_1')
zlabel('Target Function')