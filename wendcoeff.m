function [coeff, expon]=wendcoeff(d,k)
    % calculates coefficients and exponent for polynomial
    % part p_{d,k}(r) of Wendland functions
    % phi_{d,k}(r)=p_{d,k}(r)*(1-r)^expon
    expon=floor(d/2)+k+1;
    coeff= zeros(k+1,1);
    coeff(1,1)=1;
    for n=0:k-1
        coeff(n+2,1)=coeff(n+1,1)/(n+expon+2);
        for j=n+1:-1:2
            coeff(j,1)=(j*coeff(j+1,1)+coeff(j-1,1))/(expon+j);
        end
        expon=expon+1;
        coeff(1,1)=coeff(2,1)/expon;
    end
end