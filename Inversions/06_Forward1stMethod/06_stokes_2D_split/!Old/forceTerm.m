function f = forceTerm(xc) 
% 
%
f = [0; 0];

% force term for the analytic stokes problem
% x = xc(1); 
% y = xc(2);
% f = zeros(2,1);
% f(1) = 1+(-4+(12-8*y)*y)*y+(-2+(24+(-72+48*y)*y)*y+(12+(-48+(72-48*y)*y)*y+(-24+48*y+(12-24*y)*x)*x)*x)*x;
% f(2) = -12*y^2*(1-y)^2+(4+(-24+(48+(-48+24*y)*y)*y)*y+(-12+(72-72*y)*y+(8+(-48+48*y)*y)*x)*x)*x;
% f(1) = 1-2*xg-4*yg+12*yg^2-8*yg^3+24*xg*yg-72*xg*yg^2+48*xg*yg^3-48*xg^2*yg+72*xg^2*yg^2-48*xg^2*yg^3+12*xg^2-24*xg^3+48*xg^3*yg+12*xg^4-24*xg^4*yg;
% f(2) = -12*yg^2+24*yg^3-12*yg^4+48*xg*yg^2-48*xg*yg^3+24*yg^4*xg+4*xg-24*xg*yg-12*xg^2+72*xg^2*yg-72*xg^2*yg^2+8*xg^3-48*xg^3*yg+48*xg^3*yg^2; 