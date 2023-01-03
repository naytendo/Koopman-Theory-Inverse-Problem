rho1 = 1;
rho2 = 1/3;
targetType = 'matern52';
targetHyper = 1;
[interps,coefs] = getTargetFunctionExt(rho1,rho2,targetType,targetHyper);

%%


delta = [0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03];
error = zeros(1,length(delta));
condNum = zeros(1,length(delta));
for dd = 1:length(delta)
dim = 2;
% alpha0 = [89.3
% 136.75];
Gamma = 0.85;
tau = 1.25;
% 
% [boolean,kProb,checks] = isDiophantine(delta(dd),alpha0,Gamma,tau);
% 
% alpha = alpha0/norm(alpha0);



q = sqrt((0.5/delta(dd))^2-1);
alpha0 = [q
1];
% tau = 1.25;
% 
% [boolean,kProb,checks] = isDiophantine(delta,alpha0,gamma,tau);

alpha = alpha0/norm(alpha0);

boolean = 1;
if boolean
%     T = (1+dim^2*factorial(dim))^(tau+1)/Gamma/delta(dd);
    T = sqrt(q^2+1);
    t = 1:0.001:25;
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
%     plot3(x_Orb,y_Orb,z_Orb,'k');

%     target_Func_Orb = zeros(length(x_Orb),1);

%     for jj = 1:length(target_Func_Orb)
%     kernVector = zeros(length(interps),1);
%         for cc = 1:length(interps)
%             kernVector(cc) = kernel(targetType,interps(cc,:),[x_Orb(jj),y_Orb(jj),z_Orb(jj)],targetHyper);
%         end
%     target_Func_Orb(jj) = coefs'*kernVector;
%     end

    fill_Index = 1;
    while t(fill_Index) <= T && fill_Index < length(t)
        fill_Index = fill_Index+1;
    end

    xSamp = [x_Orb;y_Orb;z_Orb]';
    centers = zeros(length(xSamp),3);
    target_Func_Orb = zeros(length(xSamp),1);
    M = 1;
    centers(M,:) = xSamp(1,:);
    targetBasis = zeros(length(interps),1);
        for cc = 1:length(interps)
            targetBasis (cc) = kernel(targetType,interps(cc,:),xSamp(1,:),targetHyper);
        end
    target_Func_Orb(M) = coefs'*targetBasis ;
    separation = delta(dd);
    for mm = 1:fill_Index
        check = 0;
        for jj = 1:M
            if norm(centers(jj,:)-xSamp(mm,:)) > separation
                check = check +1;
            end
        end
        if check == M
            centers(M+1,:) = xSamp(mm,:);
            targetBasis = zeros(length(interps),1);
            for cc = 1:length(interps)
                targetBasis(cc) = kernel(targetType,interps(cc,:),xSamp(mm,:),targetHyper);
            end
            target_Func_Orb(M+1) = coefs'*targetBasis;
            M = M+1;
        end
    end
    % centers = xSamp(:,1:3:end)';
    centers = centers(1:M,:);
    target_Func_Orb = target_Func_Orb(1:M);

    beta = 1;
    K = zeros(M,M);
    type = 'wendland32';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end

    estimateCoefs = pinv(K)*target_Func_Orb;
    functionEstimate = zeros(length(t),1);
    target_Func_True = zeros(length(functionEstimate),1);

    for ii = 1:length(t)
        kernVector = zeros(length(centers),1);
        for cc = 1:length(centers)
            kernVector(cc) = kernel(type,centers(cc,:),xSamp(ii,:),beta);
        end
        functionEstimate(ii) = estimateCoefs'*kernVector;
        trueKV = zeros(length(interps),1);
        for cc = 1:length(interps)
            trueKV(cc) = kernel(targetType,interps(cc,:),xSamp(ii,:),targetHyper);
        end
        target_Func_True(ii) = coefs'*trueKV;

    end
    

    error(dd) = max(abs(functionEstimate - target_Func_True));
%     condNum(dd) = cond(K);
%     figure()
%     plot(functionEstimate)
%     hold on
%     plot(target_Func_True)
%     title(delta(dd))
end

end
%%

loglog(delta, error,'o-')
hold on
loglog(delta,100*delta.^2,'--')
% loglog(delta,condNum,'.-')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
% legend('$l_\infty$ norm','$h^2_{\Xi_n,M}$','interpreter','latex')

