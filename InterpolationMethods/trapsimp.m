

% Author: Hanieh Soleimani
% 
% This script approximates the value of the integral 
% by using the compiste Trapezoidal and composite Simpson.
% it also calculates the error of each rule with different
% subinterval m.
% 
% fun: is an inline function representing the integrand
% a and b: are the limits of integration 
% m : is the subinterval 

close all
clc
clear all

rng('default')

% create the given function 
% set the interval and m.
  fun =inline('4 + 5*x*sin(x)');
a = pi/2; b = 3*pi;
m =2;

%  fun = inline('sign(x)');  
%   a = -1; b = 1;
%   exact = fun(1) - fun(-1);
%  m = 100;

 %%%%%%%%%%%%%%%%%%%%%%%%  Composite Trapezoidal %%%%%%%%%%%%%%%%%%%%%%% 

% ezplot('4+ 5*x*sin(x)', [pi/2, 3*pi]), hold on
 
h = (b-a)/m;
trap_result=0;
for k=1:(m-1)
    x = a +h*k;
    trap_result = trap_result +feval(fun,x);
end

trap_result  = h *(feval(fun,a) +feval(fun,b))/2+ h*trap_result;
fprintf('The Result for Composite Trapezoidal:  %12.6f \n', trap_result)


% The error for the composite trapezoidal is calculated as:
%  (-1/12)*(b-a)f'')*h.^2
trap_error = (-1/12)*(b-a)*(10*cos(x) - (5*x*sin(x)))*h.^2;
fprintf('The error for the trapezoidal rule:  %12.6f \n', trap_error)

 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%  Composite Simpson %%%%%%%%%%%%%%%%%%%%%%
 
% figure, ezplot('4+ 5*x*sin(x)', [pi/2, 3*pi]), hold on
h = (b-a)/(2*m);
x1 = 0 ;
x2 = 0;
for k=1:m
    x = a + h*(2*k-1);
    x1 =x1 + feval(fun,x);
end

for k=1:(m-1)
    x= a + h*2*k;
    x2 = x2 +feval(fun,x);
end
simpson_result = h * (feval(fun,a) + feval(fun,b) + 4*x1 + 2* x2)/3;
fprintf('The Result for Composite Simpson:  % 10.6f \n', simpson_result)


%%%%%%%%%%  error for Simpson
% The error for the composite Simpson is calculated as:
% (-1/180)*(b-a)f''''] *h.^4  
simp_error = ((-1/180)*(b-a)*(5*x*sin(x) - 20*cos(x)))*h.^4;
fprintf('The error for the simpson rule:  %12.6f \n', simp_error)


 
 
 
 
 
 
 
 
 
 
 
