tspan = [0 110];
phi0 = [1;1;1];

[t,phi] = ode45(@lorenzSys,tspan,phi0);

G = @(x,y)-10*sin(1/10*y)+1/1000*(y+x).^3+200;

xSamp = phi(:,1:2);

ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = G(xSamp(ii,1),xSamp(ii,2));
end

depth = 4;
numChild = 2;
maxGrid = [30,30];
minGrid = [-30,-30];

approx = zeros(1,length(xSamp));

centers = zeros(length(xSamp),2);
M = 1;
centers(M,:) = xSamp(1,:);
separation = 0.25;
indexes = zeros(1,length(xSamp));
indexes(1) = 1;
for mm = 1:length(xSamp)
    check = 0;
    for jj = 1:M
        if norm((centers(jj,:)-xSamp(mm,:))') > separation
            check = check +1;
        end
    end
    if check == M
        centers(M+1,:) = xSamp(mm,:);
        M = M+1;
        indexes(M) = mm;
    end
end
% centers = xSamp(:,1:3:end)';
centers = centers(1:M,:);
indexes = indexes(1:M);

beta = 5;
K = zeros(M,M);
Kf = zeros(M,M);
type = 'matern32';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end

regressionCoefs = ySamp(indexes)*pinv(K);
res = 2;
x1Range = minGrid(1):res:maxGrid(1);
x2Range = minGrid(2):res:maxGrid(2);

[X1,X2] = meshgrid(x1Range,x2Range);



kernEstimate = zeros(size(X1));

for ii = 1:size(kernEstimate,1)
    for jj = 1:size(kernEstimate,2)
        kernVector = zeros(length(centers),1);
        for kk = 1:length(centers)
            kernVector(kk) = kernel(type,centers(kk,:),[X1(ii,jj),X2(ii,jj)],beta);
        end
        kernEstimate(ii,jj) = regressionCoefs*kernVector;
    end
end
%%
surf(X1,X2,G(X1,X2)- kernEstimate,'EdgeColor','none','FaceAlpha',0.4);
grid on
hold on
% surf(X1,X2,,'FaceAlpha',0.2,'FaceColor','blue','EdgeColor','none')
spacing = 4;  % play around so it fits the size of your data set
for i = 1 : spacing : length(X1(:,1))
    plot3(X1(:,i), X2(:,i), G(X1(:,i),X2(:,i))-kernEstimate(:,i),'-','color',[0.5 0.5 0.5]);
    plot3(X1(i,:), X2(i,:),G(X1(i,:),X2(i,:))- kernEstimate(i,:),'-','color',[0.5 0.5 0.5]);
end
plot3(xSamp(:,1),xSamp(:,2),zeros(length(xSamp),1),'k', 'LineWidth', 0.25)
% plot3(xSamp(:,1),xSamp(:,2),ySamp,'LineWidth', 0.25,'color','red')
% plot3(xSamp(:,1), xSamp(:,2),ySamp, 'r.', 'LineWidth', 1, 'MarkerSize', 5);

xlabel('$x$','interpreter','latex','fontsize',20)
ylabel('$y$','fontsize',20,'interpreter','latex')
zlabel('$\|(G-\hat{g}_N)(x,y)\|$','fontsize',20,'interpreter','latex')
text(12,-10,0,'$\Gamma_{\phi_0}$','fontsize',20,'interpreter','latex','color','black')
% text(-20,20,270,'$G$','fontsize',20,'interpreter','latex','color','green')
% text(-30,25,40,'$\hat{g}_N$','fontsize',20,'interpreter','latex','color','blue')
view(-24,8)
