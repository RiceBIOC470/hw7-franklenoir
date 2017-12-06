%HW7
%GB Answer
1a 70 There are two fixed points. One at 0 because a population of 0 people can’t reproduce. Any fluctuation outside of that are stably fixed at 1 because this is the carrying capacity for that population of people. At the point, the rate of growth is greatly reduced
1b. 70 you are correct, ‘a’ does not change the stability of the system. However, it does determine how quickly perturbations in X will correct towards the stable point. 
1c. 100
1d. 90 no axis labels and no discussion on how the number of fixed points evolves. 
2a. 30 This is the incorrect equation and reflects a model of gene activation and not a toggle switch. Something like this would be accepted : [V/(1+x(2)^4)-x(1); V/(1+x(1)^4)-x(2)];
2b. 95 Graphs of plotted incorrectly because your equations are wrong, but will give credit because its implemented correctly.  However, not axis labels are shown 
2c. 95 Same as 2b. 
overall: 79


%Walter Frank Lenoir 

% Problem 1: Modeling population growth
%Walter Frank Lenoir

% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).  

% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.

% Humans are unable to asexually reproduce, therefore at 0 & 1 no
% population growth will occur. 

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 

% The stability at points x = 0 and x = 1 will always be stable regardless
% of a. If x = 0, then dx/dt = a*0*(1-0/N) = 0, therefore no growth change
% is occuring. If x = 1, then dx/dt = a*1*(1-1/N) = a*(0/N) = 0. In both
% situations dx/dt = 0 at the specific points, and no growth change is occuring. 

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 

a = 4;
xo = 2;

[timecourse,time] = poptime(a,xo); %time = 1 - length(timecourse) because 
% timecourse(1) = xo. N = 1000.

% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 

clear all;

figure;
for a = 0.1:0.1:4
    for n = 1:200
    xo = rand();
    xf = xo;
    for i = 1:10 
        xf = a*xf*(1-xf);
    end
    plot(a,xf,'.');
    hold on;
    end
end

%As a approaches 2.5, xf achieves a steady state around 0.6. Once a passes
%~2.75, bifurcation patterns appear, as xf is no longer in a steady state. 

% Problem 2. Genetic toggle switches. 
% Walter Frank Lenoir

% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 

% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 

% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other

% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 

% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 

% Taken from lecture 22, dx/dt = (ku + kb(R/K)^n)/(1+(R/K)^n) - x, ku = 0, kb = V,
% R/K = A or B, x = B or A

% dA/dt = V*B^4/((1+B^4))-A, dB/dt = V*A^4/((1+A^4))-B
%
% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 
clear all;
V = 5;
rhs = @(t,x) [((V*x(2)^4)/(1+x(2)^4))-x(1); ((V*x(1)^4)/(1+x(1)^4))-x(2)];

%B0 > A0
sol = ode23(rhs, [0 10],[2,4]);
figure;
plot(sol.x, sol.y(1,:)); 
hold on;
plot(sol.x, sol.y(2,:)); 

%A0 > B0
sol = ode23(rhs, [0 10],[8,4]);
figure;
plot(sol.x, sol.y(1,:)); 
hold on;
plot(sol.x, sol.y(2,:)); 


% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter. 

gxfunc = @(x,V) ((V*x^4)/(1+x^4))-x;
for V = 1:0.1:5
    gxfunc2 = @(x) gxfunc(x,V);
    for xo = 1:0.1:5
        [rt,~,exitflag] = fzero(gxfunc2,xo);
        if exitflag == 1
            plot(V,rt,'k.');
            hold on;
        end
    end
end

