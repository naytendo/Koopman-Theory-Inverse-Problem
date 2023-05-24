function samples = samplesFromTorus(R,r)
%% Generate interpolating inputs for target function
m = 10;
samples = zeros(m,3);
theta1 = linspace(0,2*pi,m+1);
theta2 = linspace(0,2*pi,m+1);
theta1 = theta1(1:end-1);
theta2 = theta2(1:end-1);
[Theta1,Theta2] = meshgrid(theta1,theta2);
u = (R+r*cos(Theta2)).*cos(Theta1);
v = (R+r*cos(Theta2)).*sin(Theta1);
w = r*sin(Theta2);

x = zeros(m*m,3);
index = 1;
for ii = 1:m
    for jj = 1:m
        x(index,:) = [u(ii,jj),v(ii,jj),w(ii,jj)];
        index = index + 1;
    end
end
s = 2;
x = x';
for ii = 1:1000
    for jj = 2:length(x)-1
        xi = [x(:,1:jj-1),x(:,jj+1:end)];
        diff = x(:,jj)-xi;
        
        for kk = 1:length(diff)
            x(:,jj) = x(:,jj) + s*(diff(:,kk)/norm(diff(:,kk)));
        end
        
    end

end
figure()
plot3(x(1,:),x(2,:),x(3,:),'k.')


    
end



