gamma = 0.01;
gammaStr = sprintf('%g',gamma);
tend = 500;
phi0 = [1;1;1];
tspan = [0 tend];
[centers,xSamp,tNi] = generateCenters('along orbit',1,phi0,tspan);

G = @(x,y)(-10*sin(1/10*x)+1/2000*(x+y).^3)+100;
% figure()
% plot(centers(:,1),centers(:,2),'o')
% title('centers in domain')
centers = centers(1:end-1,:);
hyper = 2;
hyperStr = sprintf('%g',hyper);
M = length(centers);
K = zeros(M,M);
type = 'matern32';
    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),hyper);
        end
    end
 
%%
w0 = tNi(2:end)-tNi(1:end-1);
w = [w0(1),2*w0(2:end-1)',w0(end)]/2;
W = diag(w);
L = chol(K,'lower');
y = zeros(length(centers),1);
perturbation = zeros(length(centers),1);
delta = 5;
deltaStr = sprintf('%g',delta);

for ii = 1:length(centers)
    y(ii) = G(centers(ii,1),centers(ii,2));
    perturbation(ii) = delta*rand;
end

alpha = (K*W'*K + gamma*K)\(K*W'*y +perturbation);
beta = (L\K*W'*K/L' + gamma*eye(size(K,1),size(K,2)))\(L\K*W'*y+perturbation);

condK = cond(K*W'*K + gamma*K);
condKStr = sprintf('%g',condK);

condLK = cond(L\K*W'*K/L' + gamma*eye(size(K,1),size(K,2)));
condLKStr = sprintf('%g',condLK);

maxGrid = [30,30];
minGrid = [-30,-30];
regressionCoefsKern = alpha';
regressionCoefsNewt = beta';
res = 0.5;
x1Range = minGrid(1):res:maxGrid(1);
x2Range = minGrid(2):res:maxGrid(2);

[X1,X2] = meshgrid(x1Range,x2Range);


%%
kernEstimate = zeros(size(X1));
newtEstimate = zeros(size(X1));
for ii = 1:size(kernEstimate,1)
    for jj = 1:size(kernEstimate,2)
        kernVector = zeros(length(centers),1);
        for kk = 1:length(centers)
            kernVector(kk) = kernel(type,centers(kk,:),[X1(ii,jj),X2(ii,jj)],hyper);
        end
        kernEstimate(ii,jj) = regressionCoefsKern*kernVector;
        newtEstimate(ii,jj) = regressionCoefsNewt*(L\kernVector);
    end
end

%%
kernApprox = zeros(length(xSamp),1);
newtApprox = zeros(length(xSamp),1);
for ii = 1:length(xSamp)
        kernVector = zeros(length(centers),1);
        for kk = 1:length(centers)
            kernVector(kk) = kernel(type,centers(kk,:),xSamp(ii,:),hyper);
        end
        kernApprox(ii) = regressionCoefsKern*kernVector;
        newtApprox(ii) = regressionCoefsNewt*(L\kernVector);
    
end
%%
[kernError,maxKernIndex] = max(ySamp'-kernApprox);
[newtError,maxNewtIndex] = max(ySamp'-newtApprox);


%%
figure()
surf(X1,X2,G(X1,X2),'EdgeColor','none','FaceAlpha',0.01,'FaceColor','green');
grid on
hold on

surf(X1,X2,kernEstimate,'EdgeColor','none','FaceAlpha',0.6)
spacing2 = 3;  % play around so it fits the size of your data set

gridColor2 = [0.4 0.4 0.4];
for i = 1 : spacing2 : length(X1(:,1))
    plot3(X1(:,i), X2(:,i), kernEstimate(:,i),'-','color',gridColor2);
    plot3(X1(i,:), X2(i,:), kernEstimate(i,:),'-','color',gridColor2);
end
gridColor = [0.1 0.1 0.1];
gridSpacing = 40;  % play around so it fits the size of your data set
for i = 1 : gridSpacing : length(X1(:,1))
    plot3(X1(:,i), X2(:,i), G(X1(:,i),X2(:,i)),'-','color',gridColor);
    plot3(X1(i,:), X2(i,:), G(X1(i,:),X2(i,:)),'-','color',gridColor);
end
plot3(xSamp(:,1),xSamp(:,2),zeros(length(xSamp(:,1)),1),'k', 'LineWidth', 0.25)

ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = G(xSamp(ii,1),xSamp(ii,2));
end
plot3(xSamp(:,1), xSamp(:,2),ySamp, 'r', 'LineWidth', 0.5, 'MarkerSize', 1);


view(11,30)
ax = gca;
font_sz_ticks = 28;
font_sz_labels = 48;
zticks(200:200:400);
zlim([0,400])
ax.FontSize = font_sz_ticks;
xlabel('$x$','interpreter','latex','fontsize',font_sz_labels)
ylabel('$y$','fontsize',font_sz_labels,'interpreter','latex')
zlabel('$\hat{g}_N$','fontsize',36,'interpreter','latex')
text(10,-10,-15,'$\phi(t)$','fontsize',36,'interpreter','latex','color',[0 0 0])
% text(-30,30,450,strcat('$\gamma = $',figstr),'fontsize',font_sz_labels,'interpreter','latex')
text(0,20,300,strcat('$\gamma  =$',gammaStr),'fontsize',28,'interpreter','latex')
text(-20,20,300,strcat('$\sigma  =$',hyperStr),'fontsize',28,'interpreter','latex')
text(-20,20,370,strcat('cond # =',condKStr),'fontsize',28)
text(-40,20,300,strcat('\delta =',deltaStr),'fontsize',28)
text(20,20,300,type,'fontsize',28)
text(-20,25,200,'$G$','fontsize',font_sz_labels,'interpreter','latex','color',gridColor)
set(gcf,'Position',[100 100 800 600])
text(20,20,380,strcat('error  =',num2str(kernError)),'fontsize',28)
% plot3(xSamp(maxKernIndex,1),xSamp(maxKernIndex,2),ySamp(maxKernIndex),'ro','markerSize',20,'linewidth',4)


figure()
surf(X1,X2,G(X1,X2),'EdgeColor','none','FaceAlpha',0.01,'FaceColor','green');
grid on
hold on
surf(X1,X2,newtEstimate,'EdgeColor','none','FaceAlpha',0.6)
for i = 1 : spacing2 : length(X1(:,1))
    plot3(X1(:,i), X2(:,i), newtEstimate(:,i),'-','color',gridColor2);
    plot3(X1(i,:), X2(i,:), newtEstimate(i,:),'-','color',gridColor2);
end

gridColor = [0.1 0.1 0.1];
gridSpacing = 40;  % play around so it fits the size of your data set
for i = 1 : gridSpacing : length(X1(:,1))
    plot3(X1(:,i), X2(:,i), G(X1(:,i),X2(:,i)),'-','color',gridColor);
    plot3(X1(i,:), X2(i,:), G(X1(i,:),X2(i,:)),'-','color',gridColor);
end

plot3(xSamp(:,1),xSamp(:,2),zeros(length(xSamp(:,1)),1),'k', 'LineWidth', 0.25)

ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = G(xSamp(ii,1),xSamp(ii,2));
end
plot3(xSamp(:,1), xSamp(:,2),ySamp, 'r', 'LineWidth', 0.5, 'MarkerSize', 1);

view(11,30)
ax = gca;
font_sz_ticks = 28;
font_sz_labels = 48;
zticks(200:200:400);
zlim([0,400])
ax.FontSize = font_sz_ticks;
xlabel('$x$','interpreter','latex','fontsize',font_sz_labels)
ylabel('$y$','fontsize',font_sz_labels,'interpreter','latex')
zlabel('$\hat{g}_N$','fontsize',36,'interpreter','latex')
text(10,-10,-15,'$\phi(t)$','fontsize',36,'interpreter','latex','color',[0 0 0])
text(-20,20,300,strcat('$\sigma  =$',hyperStr),'fontsize',28,'interpreter','latex')
text(-20,20,370,strcat('cond # =',condLKStr),'fontsize',28)
text(-40,20,300,strcat('\delta =',deltaStr),'fontsize',28)
text(20,20,300,type,'fontsize',28)

% text(-30,30,450,strcat('$\gamma = $',figstr),'fontsize',font_sz_labels,'interpreter','latex')
text(0,20,300,strcat('$\gamma  =$',gammaStr),'fontsize',28,'interpreter','latex')
text(-20,25,200,'$G$','fontsize',font_sz_labels,'interpreter','latex','color',gridColor)
set(gcf,'Position',[100 100 800 600])
text(20,20,380,strcat('error  =',num2str(newtError)),'fontsize',28)
% plot3(xSamp(maxNewtIndex,1),xSamp(maxNewtIndex,2),ySamp(maxNewtIndex),'ro','markerSize',20,'linewidth',4)


