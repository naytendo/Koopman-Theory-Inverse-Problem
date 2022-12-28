rho1 = 1;
rho2 = 1/3;
targetType = 'matern52';
targetHyper = 1;
[interps,coefs] = getTargetFunction(rho1,rho2,targetType,targetHyper);

%%
delta = [0.07 0.06 0.05 0.04 0.033 0.03];
error = zeros(1,length(delta));

for dd = 1:length(delta)
dim = 2;
alpha0 = [89.3
136.75];
Gamma = 0.85;
tau = 1.25;

[boolean,kProb,checks] = isDiophantine(delta(dd),alpha0,Gamma,tau);

alpha = alpha0/norm(alpha0);

if boolean
    T = (1+dim^2*factorial(dim))^(tau+1)/Gamma/delta(dd);
    t = 1:0.1:T;
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

    xSamp = [x_Orb;y_Orb;z_Orb]';
    mod_thetaSamp = [mod_theta1_Orb,mod_theta2_Orb];
    centers = zeros(length(xSamp),2);
    target_Func_Orb = zeros(length(xSamp),1);
    M = 1;
    centers(M,:) = mod_thetaSamp(1,:);
    tfBasis = zeros(length(interps),1);
        for cc = 1:length(interps)
            tfBasis(cc) = kernel(targetType,interps(cc,:),mod_thetaSamp(1,:),targetHyper);
        end
    target_Func_Orb(M) = coefs'*tfBasis;
    separation = delta(dd);
    for mm = 1:length(mod_thetaSamp)
        check = 0;
        for jj = 1:M
            if norm(centers(jj,:)-mod_thetaSamp(mm,:)) > separation
                check = check +1;
            end
        end
        if check == M
            centers(M+1,:) = mod_thetaSamp(mm,:);
            tfBasis = zeros(length(interps),1);
            for cc = 1:length(interps)
                tfBasis(cc) = kernel(targetType,interps(cc,:),mod_thetaSamp(mm,:),targetHyper);
            end
            target_Func_Orb(M+1) = coefs'*tfBasis;
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
    functionEstimate = zeros(16000,1);
    target_Func_True = zeros(length(functionEstimate),1);

    for ii = 1:length(functionEstimate)
        estKV = zeros(length(centers),1);
        for cc = 1:length(centers)
            estKV(cc) = kernel(type,centers(cc,:),mod_thetaSamp(ii,:),beta);
        end
        functionEstimate(ii) = estimateCoefs'*estKV;
        trueBasis = zeros(length(interps),1);
        for cc = 1:length(interps)
            trueBasis(cc) = kernel(targetType,interps(cc,:),mod_thetaSamp(ii,:),targetHyper);
        end
        target_Func_True(ii) = coefs'*trueBasis;

    end
    

    error(dd) = max(abs(functionEstimate - target_Func_True));
%     figure()
%     plot(target_Func_True,'linewidth',2)
%     hold on
%     plot(functionEstimate,'--')
%     title(delta(dd))
end

end
figure()
loglog(delta, error,'o-')
hold on
loglog(delta,1500*delta.^4,'--')



