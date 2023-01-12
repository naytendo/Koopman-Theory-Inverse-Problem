rho1 = 1;
rho2 = 1/3;
targetType = 'matern';
b = 4;
dim = 3;
nu = 1/2*(b+dim/2);
targetHyper = 1;
[interps,coefs] = getTargetFunctionExt(rho1,rho2,targetType,targetHyper,nu-dim/2);

%% Creating orbit to generate centers

delta = 0.02;
scaling = sqrt((2*pi*rho1)^2+(2*pi*rho2)^2);
% scaling = 2*pi*rho2;
fillDist = 2*delta/scaling;
q = floor(sqrt((1/2/fillDist)^2-1))+1;
alpha0 = [q
1];
alpha = alpha0/norm(alpha0);
t = 1:0.01:sqrt(q^2+1)+1;
w1 = alpha(1);
w2 = alpha(2);
x = (rho1+rho2*cos(2*pi*w2*t)).*cos(2*pi*w1*t);
y = (rho1+rho2*cos(2*pi*w2*t)).*sin(2*pi*w1*t);
z = -rho2*sin(2*pi*w2*t);
xSamp = [x;y;z]'; 


%% Creating evaluation points

n1 = 256+1;
n2 = 64+1;
theta1Range = linspace(0,2*pi,n1);
theta2Range = linspace(0,2*pi,n2);
theta1Range = theta1Range(1:end-1);
theta2Range = theta2Range(1:end-1);
Eval_points = (n1-1)*(n2-1);
[theta1Evals,theta2Evals] = meshgrid(theta1Range,theta2Range);


u = (rho1+rho2*cos(theta2Evals)).*cos(theta1Evals);
v = (rho1+rho2*cos(theta2Evals)).*sin(theta1Evals);
w = rho2*sin(theta2Evals);

x_Eval = zeros(Eval_points,3);
index = 1;
for ii = 1:n2-1
    for jj = 1:n1-1
        x_Eval(index,:) = [u(ii,jj),v(ii,jj),w(ii,jj)];
        index = index + 1;
    end
end

%%


% delta = [0.045 0.038 0.033 0.025];
delta = [0.4 0.2 0.15 0.12 0.1 0.08 0.06 0.04 0.03 0.025 0.02];
linfty_error = zeros(1,length(delta));
l2_error = zeros(1,length(delta));
condNum = zeros(1,length(delta));

for dd = 1:length(delta)


boolean = 1;
alpha = alpha0/norm(alpha0);
if boolean


    centers = zeros(length(xSamp),3);
    target_Func_Orb = zeros(length(xSamp),1);
    M = 1;
    centers(M,:) = xSamp(1,:);
    targetBasis = zeros(size(interps,1),1);
        for cc = 1:size(interps,1)
            targetBasis (cc) = kernel(targetType,interps(cc,:),xSamp(1,:),targetHyper,nu-dim/2);
        end
    target_Func_Orb(M) = coefs'*targetBasis ;
    separation = delta(dd);
    for mm = 1:length(xSamp)
        check = 0;
        for jj = 1:M
            if norm(centers(jj,:)-xSamp(mm,:)) > 2*separation
                check = check +1;
            end
        end
        if check == M
            centers(M+1,:) = xSamp(mm,:);
            targetBasis = zeros(length(interps),1);
            for cc = 1:length(interps)
                targetBasis(cc) = kernel(targetType,interps(cc,:),xSamp(mm,:),targetHyper,nu-dim/2);
            end
            target_Func_Orb(M+1) = coefs'*targetBasis;
            M = M+1;
        end
    end
    % centers = xSamp(:,1:3:end)';
    centers = centers(1:M,:);
    target_Func_Orb = target_Func_Orb(1:M);

    hyperParam = 8/3;
    K = zeros(M,M);
    type = 'wendland32';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),hyperParam,[]);
        end
    end
%     L = chol(K,'lower');
    estimateCoefs = (K)\(target_Func_Orb);
    

    
    functionEstimate = zeros(Eval_points,1);
    target_Func_True = zeros(length(functionEstimate),1);


%%
    for ii = 1:Eval_points
        kernVector = zeros(size(centers,1),1);
        for cc = 1:size(centers,1)
            kernVector(cc) = kernel(type,centers(cc,:),x_Eval(ii,:),hyperParam,[]);
        end
        functionEstimate(ii) = estimateCoefs'*kernVector;
        trueKV = zeros(length(interps),1);
        for cc = 1:length(interps)
            trueKV(cc) = kernel(targetType,interps(cc,:),x_Eval(ii,:),targetHyper,nu-dim/2);
        end
        target_Func_True(ii) = coefs'*trueKV;

    end
    error = abs(functionEstimate - target_Func_True); 
    

    linfty_error(dd) = max(error);
    l2_error(dd) = sqrt(sum(error.^2));
    condNum(dd) = cond(K);
%     figure()
%     plot(functionEstimate)
%     hold on
%     plot(target_Func_True)
%     title(delta(dd))
% figure()
% plot(mod_theta2_Orb,mod_theta1_Orb,'.')
% title(strcat('\delta =',num2str(delta(dd)),', Direction = [',num2str(q), ',1]^T','.'))

figure()
plot3(centers(:,1),centers(:,2),centers(:,3),'ko')
hold on
plot3(x,y,z,'k')
title(strcat('\delta =',num2str(delta(dd))))
end

end
%%
figure()
loglog(delta, linfty_error,'rx-')
hold on
loglog(delta, 0.01*l2_error,'bo-')
loglog(delta,8^2*delta.^(2.5),'k-.')
loglog(delta,5^3*delta.^(3.5),'k--')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
ylabel('Error')
 legend('$l_\infty$ error','$l_2$ error','$h^{2.5}_{\Xi_n,M}$','$h^{3.5}_{\Xi_n,M}$','interpreter','latex')
figure()
loglog(delta,condNum,'.-')
ylabel('cond(K)')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')


