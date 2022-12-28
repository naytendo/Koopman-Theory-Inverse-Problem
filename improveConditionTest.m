tspan = [0 110];
phi0 = [1;1;1];

[t,phi] = ode45(@lorenzSys,tspan,phi0);

% G = @(x,y)-10*sin(1/10*y)+1/1000*(y+x).^3+200;
sep = [0.5,0.1,0.05,0.01];
condK = zeros(1,length(sep));
condLKL = zeros(1,length(sep));
minEigK = zeros(1,length(sep));
minEigLKL = zeros(1,length(sep));
for ii = 1:length(sep)
[centers,xSamp,tNi] = generateCenters('along orbit',sep(ii),phi0,tspan);

centers = centers(1:end-1,:);
beta = 0.25;
M = length(centers);
K = zeros(M,M);
type = 'matern32';
    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end

L = chol(K,'lower');

w0 = tNi(2:end)-tNi(1:end-1);
w = [w0(1),2*w0(2:end-1)',w0(end)]/2;
W = diag(w);

gamma = 1;
condK(ii) = cond(K*W'*K + gamma*eye(size(K,1),size(K,2)));
condLKL(ii) = cond(L\K*W'*K/L' + gamma*eye(size(K,1),size(K,2)));
minEigK(ii) = min(eig(K*W'*K + gamma*eye(size(K,1),size(K,2))));
minEigLKL(ii) = min(eig(L\K*W'*K/L' + gamma*eye(size(K,1),size(K,2))));
end

%%
plot(-log(sep),log(condK),'-o')

hold on
plot(-log(sep),log(condLKL),'-.')
ylabel('log(condition #)')
xlabel('-log(sep)')
plot(-log(sep),log(1./minEigK),'-s')
plot(-log(sep),log(1./minEigLKL),'-+')