rho1 = 1;
rho2 = 1/3;
targetType = 'matern';
b = 4;
dim = 3;
nu = 1/2*(b+dim/2);
targetHyper = 1;
[interps,coefs] = getTargetFunctionExt(rho1,rho2,targetType,targetHyper,nu-dim/2);

%% Creating orbit to generate centers

delta = 0.01;
scaling = sqrt((2*pi*rho1)^2+(2*pi*rho2)^2);
% scaling = 2*pi*rho2;
fillDist = 2*delta/scaling;
q = floor(sqrt((1/2/fillDist)^2-1))+1;
alpha0 = [q
1];
alpha = alpha0/norm(alpha0);
t = 0:0.01:sqrt(q^2+1);
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
mu = rho2*sin(theta2Evals);

x_Eval = zeros(Eval_points,3);
index = 1;
for ii = 1:n2-1
    for jj = 1:n1-1
        x_Eval(index,:) = [u(ii,jj),v(ii,jj),mu(ii,jj)];
        index = index + 1;
    end
end

%% Determine Quadrature Points
h = 0.05;
Fints = zeros(length(xSamp),3);
for mm = 1:length(xSamp)
    check = 0;
    for jj = 1:M
        if norm(centers(jj,:)-xSamp(mm,:)) > 2*separation
            check = check +1;
        end
    end
    if check == M
        Fints(M+1,:) = xSamp(mm+1,:);
        M = M+1;
    end
end
Fints = Fints(1:M,:);

mu0 = zeros(size(Fints,1)-1,1);
fquads = zeros(size(Fints,1)-1,size(Fints,2));
for ii = 1:length(Fints)-1
    mu0(ii) = norm(Fints(ii+1,:)-Fints(ii,:));
    fquads(ii,:) = (Fints(ii+1,:)+Fints(ii,:))/2;
end
mu = [mu0(1);2*mu0(2:end-1);mu0(end)]/2;
W = diag(mu);



%%

% delta = [0.045 0.038 0.033 0.025];
delta = [0.4 0.2 0.15 0.12 0.1 0.08 0.06 0.04];
linfty_error = zeros(1,length(delta));
l2_error = zeros(1,length(delta));
condNum = zeros(1,length(delta));
condNumRed = zeros(1,length(delta));
condNumReg = zeros(1,length(delta));




for dd = 1:length(delta)


alpha = alpha0/norm(alpha0);



centers = zeros(length(xSamp),3);

target_Func_Orb = zeros(length(xSamp),1);
N = 1;
centers(N,:) = xSamp(1,:);
targetBasis = zeros(size(interps,1),1);
    for cc = 1:size(interps,1)
        targetBasis (cc) = kernel(targetType,interps(cc,:),xSamp(1,:),targetHyper,nu-dim/2);
    end
target_Func_Orb(N) = coefs'*targetBasis ;
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
        targetBasis = zeros(length(interps),1);
        for cc = 1:length(interps)
            targetBasis(cc) = kernel(targetType,interps(cc,:),xSamp(nn,:),targetHyper,nu-dim/2);
        end
        target_Func_Orb(M+1) = coefs'*targetBasis;
        N = N+1;
    end
end

    %%
% centers = xSamp(:,1:3:end)';
centers = centers(1:N,:);


target_Func_Orb = target_Func_Orb(1:N);

hyperParam = 8/3;
K = zeros(M-1,M-1);
type = 'wendland32';

for pp = 1:M
    for jj = 1:N
        K(pp,jj) = kernel(type,centers(jj,:),Fints(pp,:),hyperParam,[]);
    end
end
%     L = chol(K,'lower');
%     A = L\K;
    A = K*W*K';
    [U,S,V] = svd(A);
    sings = diag(S);
    singsR = sings(log(sings)>-5);
    r = length(singsR);
    figure()
    plot(log(sings))
    title(sprintf('SVD of A = KWK^T with fill = %.3f',delta(dd)));
    ylabel('\sigma_A')

    gamma = 0.01;
    condNum(dd) = cond(A);
    condNumRed(dd) = cond(U(:,1:r)*S(1:r,1:r)*V(:,1:r)');
    condNumReg(dd) = cond(A+gamma*eye(M-1,M-1));

end
%%
figure()
loglog(delta,condNum,'-o')
hold on
loglog(delta,condNumRed,'-x')
loglog(delta,condNumReg,'-s')

legend('Standard','Reduced SVD','Regularized')

