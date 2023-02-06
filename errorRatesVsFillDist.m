rho1 = 1/4;
rho2 = 1/6;
targetType = 'matern';
b = 4;
dim = 3;
nu = 1/2*(b+dim/2);
targetHyper = 1;
[interps,coefs] = getTargetFunctionExt(rho1,rho2,targetType,targetHyper,nu-dim/2);

%% Creating filling orbit 

delta = 0.007;
scaling = sqrt((2*pi*rho1)^2+(2*pi*rho2)^2)/sqrt(2);
% scaling = 2*pi*rho2;
fillDist = 2*delta/scaling;
q = floor(sqrt((1/2/fillDist)^2-1))+1;
alpha0 = [q
1];
alpha = alpha0/norm(alpha0);
t = 0:0.0005:sqrt(q^2+1);
w1 = alpha(1);
w2 = alpha(2);
x = (rho1+rho2*cos(2*pi*w2*t)).*cos(2*pi*w1*t);
y = (rho1+rho2*cos(2*pi*w2*t)).*sin(2*pi*w1*t);
z = -rho2*sin(2*pi*w2*t);
xSamp = [x;y;z]'; 

x_Eval = xSamp;
Eval_points = length(x_Eval);
% %% Creating evaluation points
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

%% Determine Quadrature Points

M0 = 1;
fints = zeros(size(xSamp,1)-1,size(xSamp,2));
fints(M0,:) = xSamp(2,:);
last = fints(M0,:);

h = 0.01;
for mm = 1:length(xSamp)-1
    current = xSamp(mm+1,:);
    if norm(last-current) > h
        M0 = M0+1;
        fints(M0,:) = xSamp(mm+1,:);
        last = current;
    end
end

fints = fints(1:M0,:);


M = length(fints)-1;
yQuads= zeros(M,1);
mu0 = zeros(M,1);
%%
fquads = zeros(M,size(fints,2));
for ii = 1:M
    mu0(ii) = norm(fints(ii+1,:)-fints(ii,:));
    fquads(ii,:) = (fints(ii+1,:)+fints(ii,:))/2; % taving an average
    targetBasis = zeros(length(interps),1);
    for cc = 1:length(interps)
        targetBasis(cc) = kernel(targetType,interps(cc,:),fquads(ii,:),targetHyper,nu-dim/2);
     end
     yQuads(ii) = coefs'*targetBasis;
end


W = diag(mu0);


%%

% delta = [0.045 0.038 0.033 0.025];
delta = [0.1 0.08 0.07 0.06 0.05 0.04 0.033 0.028 0.025 0.022 0.02 0.018];
linfty_error = zeros(1,length(delta));
l2_error = zeros(1,length(delta));
condNum = zeros(1,length(delta));
condNumRed = zeros(1,length(delta));
condNumReg = zeros(1,length(delta));
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
    hyperParam = 1;
    K = zeros(N,M);
    type = 'wendland32';
    for pp = 1:M
        for jj = 1:N
            K(jj,pp) = kernel(type,centers(jj,:),fquads(pp,:),hyperParam,[]);
        end
    end
    %     L = chol(K,'lower');
    %     A = L\K;
        A = K*W*K';
        estimateCoefs = (A)\(K*W*yQuads);
        functionEstimate = zeros(Eval_points,1);
        target_Func_True = zeros(length(functionEstimate),1);
    
    
        for ii = 1:Eval_points-1
            kernVector = zeros(size(centers,1),1);
            for cc = 1:size(centers,1)
                kernVector(cc) = kernel(type,centers(cc,:),x_Eval(ii+1,:),hyperParam,[]);
            end
            functionEstimate(ii) = estimateCoefs'*(kernVector);
            trueKV = zeros(length(interps),1);
            for cc = 1:length(interps)
                trueKV(cc) = kernel(targetType,interps(cc,:),x_Eval(ii+1,:),targetHyper,nu-dim/2);
            end
            target_Func_True(ii) = coefs'*trueKV;

        end
    error = abs(functionEstimate - target_Func_True); 
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

%%
figure()
loglog(delta, linfty_error,'rx-')
hold on
loglog(delta, 0.01*l2_error,'bo-')
loglog(delta,0.25*delta.^(2.5),'r-.')
loglog(delta,1.5*delta.^(3.5),'b--')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
ylabel('Error')
legend('$l_\infty$ error','$l_2$ error','$h^{2.5}_{\Xi_n,M}$','$h^{3.5}_{\Xi_n,M}$','interpreter','latex')
set(gca,'fontsize',20)
% legend('$l_\infty$ error','$l_2$ error','interpreter','latex')

%%
figure()
loglog(delta,condNum,'.-')
ylabel('cond(K)')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
set(gca,'fontsize',20)

