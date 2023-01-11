function kern = kernel(type,x,y,beta,nu)
% Calculate value of kernel given the type of kernel k(.,.),two elements
% and parameters.
% Input:
%       type    -   string, type of kernel 
%       x,y     -   column vectors
%       beta    -   column vectors, parameters of kernel (indications to be added)
% Output:
%       KK      -   k(x,y)
r = norm(x-y);
% Examine the inputs (to be added)
switch type
    case 'exp'
        Dxy = r^2;
        kern = 1/(2*pi*beta^length(x))*exp(-1/2/beta^2 * Dxy);
    case 'wendland31'
        if 1-r <= 0
            kern = 0;
        else
            kern = (1-r)^4*(4*r/beta+1);
        end
     case 'wendland32'
        if 1-r/beta <= 0
            kern = 0;
        else
            kern = (1-r/beta)^6*(35*(r/beta)^2 + 18*r/beta +3);
        end
     case 'wendland33'
        if 1-r <= 0
            kern = 0;
        else
            kern = (1-r)^8*(32*r^3/(beta^3) + 25*r^2/(beta^2) + 8*r/beta + 1);
        end
    case 'imq'
        kern = 1/(5^2+r^2)^beta;
    case 'matern32'
        kern = (1 + sqrt(3)*r/beta)*exp(-sqrt(3)*r/beta);
    case 'matern52'
        kern = (1 + sqrt(5)*r/beta + 5/3*r^2/(beta^2))*exp(-sqrt(5)*r/beta);
    case 'matern72'
        kern = (1+ sqrt(7)*r + 7/5*r^2/(beta^2)+7/15*sqrt(7)*r^3/(beta^3))*exp(-sqrt(7)*r/beta);
    case 'matern'
        if r > 1e-15
            sigma = 1;
            kern = sigma^2*2^(1-nu)/gamma(nu)*(sqrt(2*nu)*r/beta)^nu*besselk(nu,sqrt(2*nu)*r/beta);
        else
            kern = 1;
        end
    otherwise 
   
        error('Kernel type is not supported!');
end