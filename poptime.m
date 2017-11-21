function [timecourse,time] = poptime(a,xo)

if (xo > 1)

N = 1000*0.99;

rhs = @(t,x) [a*x(1)*(1-x(1))];

x = xo;
X = x/N;

sol = ode23(rhs, [x N],X);

timecourse = (sol.x);
time = length(timecourse)-1;
  
else
    disp('xo < 2')
end
end 