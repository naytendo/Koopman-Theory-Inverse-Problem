rho1 = 3;
rho2 = 1/3;
targetType = 'ms';
b = 4;
dim = 3;
nu = 1/2*(b+dim/2);
targetScale = 1;
[interps,coefs] = getTargetFunctionAlongOrbit2(rho1,rho2,targetType,targetScale,nu-dim/2);

%% Creating filling orbit 

% delta = 0.007;
% scaling = sqrt((2*pi*rho1)^2+(2*pi*rho2)^2)/sqrt(2);
% % scaling = 2*pi*rho2;
% fillDist = 2*delta/scaling;
% q = floor(sqrt((1/2/fillDist)^2-1))+1;
% alpha0 = [q
% 1];
% alpha = alpha0/norm(alpha0);
% t = 0:0.0005:sqrt(q^2+1);
% w1 = alpha(1);
% w2 = alpha(2);
% x = (rho1+rho2*cos(2*pi*w2*t)).*cos(2*pi*w1*t);
% y = (rho1+rho2*cos(2*pi*w2*t)).*sin(2*pi*w1*t);
% z = -rho2*sin(2*pi*w2*t);
% xSamp = [x;y;z]'; 

q = 1/6;
alpha0 = [q
1];
alpha = alpha0/norm(alpha0);
t = 0:0.001:sqrt((1/q)^2+1);
w1 = alpha(1);
w2 = alpha(2);
x = (rho1+rho2*cos(2*pi*w2*t)).*cos(2*pi*w1*t);
y = (rho1+rho2*cos(2*pi*w2*t)).*sin(2*pi*w1*t);
z = -rho2*sin(2*pi*w2*t);
xSamp = [x;y;z]'; 

% Creates Eval Points
x_Eval = zeros(length(xSamp),3); 
EvPts = 1;
x_Eval(EvPts,:) = xSamp(1,:);
separation = 0.025;
for nn = 1:length(xSamp)
    check = 0;
    for jj = 1:EvPts
        if norm(x_Eval(jj,:)-xSamp(nn,:)) > 2*separation
            check = check +1;
        end
    end
    if check == EvPts
        x_Eval(EvPts+1,:) = xSamp(nn,:);
        EvPts = EvPts+1;
    end
end


Eval_points = EvPts;
x_Eval = x_Eval(1:EvPts,:);
% Another way to creating evaluation points
% 
% n1 = 256+1;
% n2 = 64+1;
% theta1Range = linspace(0,2*pi,n1);
% theta2Range = linspace(0,2*pi,n2);
% theta1Range = theta1Range(1:end-1);
% theta2Range = theta2Range(1:end-1);
% Eval_points = (n1-1)*(n2-1);
% [theta1Evals,theta2Evals] = meshgrid(theta1Range,theta2Range);
% 
% 
% u = (rho1+rho2*cos(theta2Evals)).*cos(theta1Evals);
% v = (rho1+rho2*cos(theta2Evals)).*sin(theta1Evals);
% mu = rho2*sin(theta2Evals);
% 
% x_Eval = zeros(Eval_points,3);
% index = 1;
% for ii = 1:n2-1
%     for jj = 1:n1-1
%         x_Eval(index,:) = [u(ii,jj),v(ii,jj),mu(ii,jj)];
%         index = index + 1;
%     end
% end
figure(1)
hold on

figure(2)
hold on

figure(3)
hold on

mkr = ['-s','-o','-d','-v'];
%% Determine Quadrature Points
hQ = [0.1 0.08 0.04 0.01];
for hh = 1:length(hQ)
M0 = 1;
fints = zeros(size(xSamp,1)-1,size(xSamp,2));
fints(M0,:) = xSamp(2,:);
last = fints(M0,:);

for mm = 1:length(xSamp)-1
    current = xSamp(mm+1,:);
    if norm(last-current) > hQ(hh)
        M0 = M0+1;
        fints(M0,:) = xSamp(mm+1,:);
        last = current;
    end
end

fints = fints(1:M0,:);


M = length(fints)-1;
mu0 = zeros(M,1);
%%
fquads = zeros(M,size(fints,2));
for ii = 1:M
    mu0(ii) = norm(fints(ii+1,:)-fints(ii,:));
    fquads(ii,:) = (fints(ii+1,:)+fints(ii,:))/2; % taking an average
end

Kquads = kermat(interps,fquads,targetType,nu-dim/2,targetScale);
yQuads = coefs'*Kquads;
W = diag(mu0);

%%

% delta = [0.045 0.038 0.033 0.025];
delta = [0.8 0.4 0.2 0.1 0.05 0.025];
linfty_error = zeros(1,length(delta));
l2_error = zeros(1,length(delta));
condNum = zeros(1,length(delta));
% condNumRed = zeros(1,length(delta));
% condNumReg = zeros(1,length(delta));
funcDim = zeros(1,length(delta));

for dd = 1:length(delta)
    centers = zeros(length(xSamp),3); 
    N = 1;
    centers(N,:) = xSamp(1,:);
    separation = delta(dd);
    for nn = 1:length(xSamp)
        check = 0;
        for jj = 1:N
            if norm(centers(jj,:)-xSamp(nn,:)) > 2*separation
                check = check +1;
            end
        end
        if check == N
            centers(N+1,:) = xSamp(nn,:);
            N = N+1;
        end
    end
    % centers = xSamp(:,1:3:end)';
    centers = centers(1:N,:);
    scale = 1.7;
    NatPar = 2;
    type = 'w3';
%     for pp = 1:M
%         for jj = 1:N
%             K(jj,pp) = kernel(type,centers(jj,:),fquads(pp,:),hyperParam,[]);
%         end
%     end

    K = kermat(centers,fquads,type,NatPar,scale);
    %     L = chol(K,'lower');
    %     A = L\K;
    A = K*K';
    estCoefs = (A)\(K*yQuads');

    Kest = kermat(centers,x_Eval,type,NatPar,scale);

    KAct = kermat(interps,x_Eval,targetType,nu-dim/2,targetScale);
    yEst = estCoefs'*Kest;
    yAct = coefs'*KAct;
    %%
    error = abs(yEst - yAct); 
    linfty_error(dd) = max(error);
    l2_error(dd) = sqrt(sum(error.^2));
    funcDim(dd) = N;
    condNum(dd) = cond(A);
%     figure()
%     plot(functionEstimate)
%     hold on
%     plot(target_Func_True)
%     title(delta(dd))
% figure()
% plot(mod_theta2_Orb,mod_theta1_Orb,'.')
% title(strcat('\delta =',num2str(delta(dd)),', Direction = [',num2str(q), ',1]^T','.'))

% figure()
% plot3(centers(:,1),centers(:,2),centers(:,3),'ko')
% hold on
% plot3(x,y,z,'k')
% title(strcat('\delta =',num2str(delta(dd))))
end
figure(1)
loglog(delta, 0.03*l2_error,mkr(2*hh-1:2*hh),'linewidth',1.5)


figure(2)
loglog(delta, condNum,mkr(2*hh-1:2*hh),'linewidth',1.5)


figure(3)
loglog(delta, funcDim,mkr(2*hh-1:2*hh),'linewidth',1.5)


% legend('$l_\infty$ error','$l_2$ error','interpreter','latex')
end
figure(1)
loglog(delta,100*delta.^(3),'k-.')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
ylabel('$l_\infty$ error','interpreter','latex')
legend(sprintf('$h_Q = %.2f$',hQ(1)),sprintf('$h_Q = %.2f$',hQ(2)),sprintf('$h_Q = %.2f$',hQ(3)),sprintf('$h_Q = %.2f$',hQ(4)),'$h^{3}_{\Xi_n,M}$','interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')

figure(2)
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
ylabel('Condition Number','interpreter','latex')
legend(sprintf('$h_Q = %.2f$',hQ(1)),sprintf('$h_Q = %.2f$',hQ(2)),sprintf('$h_Q = %.2f$',hQ(3)),sprintf('$h_Q = %.2f$',hQ(4)),'interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')

figure(3)
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
ylabel('$N$','interpreter','latex')
legend(sprintf('$h_Q = %.2f$',hQ(1)),sprintf('$h_Q = %.2f$',hQ(2)),sprintf('$h_Q = %.2f$',hQ(3)),sprintf('$h_Q = %.2f$',hQ(4)),'interpreter','latex')
set(gca,'fontsize',20)
set(gca,'XScale','log')
set(gca,'YScale','log')