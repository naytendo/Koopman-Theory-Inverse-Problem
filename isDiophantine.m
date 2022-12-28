function [boolean,kProb,checks] = isDiophantine(delta,alpha,gamma,tau)
N = floor(9/delta)+1;
k = [-N;0];
checks = 0;
kProb = [];
boolean = 1;
% figure()
while k(1) <= N
    k(2) = -floor(sqrt(N^2-k(1)^2));
    while norm(k) <= N
        if k(1) == 0 && k(2) == 0
        elseif abs(dot(alpha,k)) < gamma*norm(k)^-tau
           kProb = k;
           boolean = 0;
           break
        end
        k(2) = k(2) + 1;
        checks = checks +1;
%         plot(k(1),k(2),'kx')
%         hold on
        
    end
    if ~isempty(kProb)
        break
    end
    k(1) = k(1) + 1;
end
   
end
        