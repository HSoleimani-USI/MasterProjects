

% Author: Hanieh Soleimani
% This script coputes the interpolation by
% using divide difference newton interpolation.
close all
clc
clear all

rng('default')


%  the degree of polynomial
n =2;
x = linspace(-1,1,n+1);
fun = 1./(2 + 3*(x.*x));
 
z = length(x);
D = zeros(z,z);
D(:,1) = fun';
% creating divided difference coefficients
for i=2:z
  for j=i:z
      D(j,i) = (D(j,i-1)-D(j-1,i-1))/(x(j)-x(j-i+1));
  end
end
C = D(z,z);
for k=(z-1):-1:1
% using matlab function: C = conv(A, B) convolves vectors A and B.
  C = conv(C,poly(x(k)));
  c = length(C);
  C(c) = C(c) + D(k,k);
end


xeval = -1:0.1:1;
yeval = polyval(C,xeval);
plot(xeval,yeval); hold on;















