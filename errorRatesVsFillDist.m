rho1 = 1;
rho2 = 1/3;
targetType = 'matern';
b = 4;
dim = 2;
nu = 0.5*(b-dim/2);
targetHyper = 1;
[interps,coefs] = getTargetFunctionExt(rho1,rho2,targetType,targetHyper,nu-dim/2);

%% Creating evaluation orbit

q_Eval = 105;
alpha0_Eval = [q_Eval,1];
alpha_Eval = alpha0_Eval/norm(alpha0_Eval);
t_Eval = 1:0.005:sqrt(q_Eval^2+1)+1;
w1_Eval = alpha_Eval(1);
w2_Eval = alpha_Eval(2);
x_E = (rho1+rho2*cos(2*pi*w2_Eval*t_Eval)).*cos(2*pi*w1_Eval*t_Eval);
y_E = (rho1+rho2*cos(2*pi*w2_Eval*t_Eval)).*sin(2*pi*w1_Eval*t_Eval);
z_E = -rho2*sin(2*pi*w2_Eval*t_Eval);
x_Eval = [x_E;y_E;z_E]'; 

%%


% delta = [0.045 0.038 0.033 0.025];
delta = [0.4 0.2 0.15 0.12 0.1 0.08 0.06 0.04 0.03 0.025 0.02];
error = zeros(1,length(delta));
condNum = zeros(1,length(delta));

for dd = 1:length(delta)

% alpha0 = [89.3
% 136.75];
% Gamma = 0.85;
% tau = 1.25;
% 
scaling = sqrt((2*pi*rho1)^2 +(2*pi*rho2)^2);
deltaFill = delta(dd)/(2*pi*rho2);
q = floor(sqrt((1/2/deltaFill)^2-1))+1;
alpha0 = [q
1];
% [boolean,kProb,checks] = isDiophantine(delta(dd),alpha0,Gamma,tau);
boolean = 1;
alpha = alpha0/norm(alpha0);
if boolean

%     T = (1+dim^2*factorial(dim))^(tau+1)/Gamma/(delta(dd)^tau);
    T = sqrt(q^2+1)+1;
    t = 1:0.01:60;
    
    
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

    theta1_Orb = alpha(1)*t(1:fill_Index);
    theta2_Orb = alpha(2)*t(1:fill_Index);
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

    xSamp = [x_Orb;y_Orb;z_Orb]';
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
    for mm = 1:fill_Index
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

    estimateCoefs = pinv(K)*target_Func_Orb;
    

    
    functionEstimate = zeros(length(t_Eval),1);
    target_Func_True = zeros(length(functionEstimate),1);


%%
    for ii = 1:length(t_Eval)
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
    

    error(dd) = max(abs(functionEstimate - target_Func_True));
    condNum(dd) = cond(K);
%     figure()
%     plot(functionEstimate)
%     hold on
%     plot(target_Func_True)
%     title(delta(dd))
figure()
plot(mod_theta2_Orb,mod_theta1_Orb,'.')
title(strcat('\delta =',num2str(delta(dd)),', Direction = [',num2str(q), ',1]^T','.'))

figure()
plot3(centers(:,1),centers(:,2),centers(:,3),'k.','markerSize',2)
title(strcat('\delta =',num2str(delta(dd)),', Direction = [',num2str(q), ',1]^T'))
end

end
%%
figure()
loglog(delta, error,'o-')
hold on
loglog(delta,10^3*delta.^(2.5),'--')

figure()
loglog(delta,condNum,'.-')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
% legend('$l_\infty$ norm','$h^2_{\Xi_n,M}$','interpreter','latex')

