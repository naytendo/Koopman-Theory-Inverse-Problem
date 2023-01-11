rho1 = 1;
rho2 = 1/3;
targetType = 'matern';
b = 4;
dim = 3;
nu = 1/2*(b+dim/2);
targetHyper = 1;
[interps,coefs] = getTargetFunctionExt(rho1,rho2,targetType,targetHyper,nu-dim/2);

%% Creating evaluation orbit

q_Eval = 120;
alpha0_Eval = [q_Eval,1];
alpha_Eval = alpha0_Eval/norm(alpha0_Eval);
t_Eval = 1:0.005:sqrt(q_Eval^2+1)+1;
w1_Eval = alpha_Eval(1);
w2_Eval = alpha_Eval(2);
x_E = (rho1+rho2*cos(2*pi*w2_Eval*t_Eval)).*cos(2*pi*w1_Eval*t_Eval);
y_E = (rho1+rho2*cos(2*pi*w2_Eval*t_Eval)).*sin(2*pi*w1_Eval*t_Eval);
z_E = -rho2*sin(2*pi*w2_Eval*t_Eval);
x_Eval = [x_E;y_E;z_E]'; 

%% Creating orbit to generate centers

delta = 0.02;
% scaling = sqrt((2*pi*rho1)^2+(2*pi*rho2)^2);
scaling = 2*pi*rho2;
deltaScaled = delta/scaling;
q = floor(sqrt((1/2/deltaScaled)^2-1))+1;
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

%%


% delta = [0.045 0.038 0.033 0.025];
delta = [0.4 0.2 0.15 0.12 0.1 0.08 0.06 0.04 0.03 0.025 0.02];
error = zeros(1,length(delta));
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
% figure()
% plot(mod_theta2_Orb,mod_theta1_Orb,'.')
% title(strcat('\delta =',num2str(delta(dd)),', Direction = [',num2str(q), ',1]^T','.'))

% figure()
% plot3(centers(:,1),centers(:,2),centers(:,3),'ko')
% hold on
% plot3(x,y,z,'k')
% title(strcat('\delta =',num2str(delta(dd)),', Direction = [',num2str(q), ',1]^T'))
end

end
%%
figure()
loglog(delta, error,'o-')
hold on
loglog(delta,10^3*delta.^(2.5),'--')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
ylabel('Error')
figure()
loglog(delta,condNum,'.-')
ylabel('cond(K)')
xlabel('$h_{\Xi_n,M}$','interpreter','latex')
% legend('$l_\infty$ norm','$h^2_{\Xi_n,M}$','interpreter','latex')

