
% Author: Hanieh Soleimani


%  Application of Euler methods,
% with both implicit and explicit one, then
% test them with the numerical approximations in a plot.

close all
clc
clear all

rng('default');


% set the value to lambda for 2 and -10
% set different step size h as given
 l = -10;
% l = 2; 
h = [.1, .05, .005]; hIndex = 1;


gridpoints = [(1-0)/h(1), (1-0)/h(2), (1-0)/h(3)];
pointIndex = 1;

%  set t with the given value between 0 and 1
%  setting the analytical one with 
%  soution to the given ODE
t = linspace(0,1,gridpoints(1));
exact = 2*exp(l*t);


%%%%%%%%%%%%%%%%%%%%%% Applying the EXPLICIT EULER method  %%%%%%%%%%%%%%
figure(1);

%  plot the exact solution in black line
plot(t, exact, 'k');

xlabel('t');
ylabel('f(t)');
title('EXACT AND APPROXIMATED EXPLICIT EULER');
hold on;

tic;
for k = 1:3
    % set the initial condition 
    explicitEuler(1,1) = 2;
    
    for i = 2:(gridpoints(pointIndex))
        explicitEuler(i,1) = (explicitEuler(i-1,1) + h(k)*l*explicitEuler(i-1,1));
    end
    
    plot(linspace(0,1,gridpoints(pointIndex)), explicitEuler, '.-');
    
    pointIndex = pointIndex+1;
    
end
toc;

legend('exact','h= .1','h= .05','h= .005','Location','southoutside','Orientation','horizontal');



%%%%%%%%%%%%%%%%%%%%%% Applying the IMPLICIT EULER method  %%%%%%%%%%%%%%

figure(2);

%  plot the exact solution in black line
plot(t, exact, 'k');

xlabel('t');
ylabel('f(t)');
title('EXACT AND APPROXIMATED IMPLICIT EULER');

hold on;

tic;
for k = 1:3
    %set the initial condition 
    implicitEuler(1,1) = 2;
    
    for i = 2:(gridpoints(hIndex))
        implicitEuler(i,1) = (implicitEuler(i-1,1) / (1 - h(k)*l));
    end
    
    plot(linspace(0,1,gridpoints(hIndex)), implicitEuler, '.-');
    
    hIndex = hIndex+1;
    
end
toc;



legend('exact','h= .1','h= .05','h= .005','Location','southoutside','Orientation','horizontal');





