% Backward Euler Methods 
clc;
clear all;
close all;
% Analysis along the time [a,b] 
a = 0; % The initial time 
b = 2; % The final time 

% The number of the "jumps" within [a,b] 
N = 30;

% The step size 
h = (b-a)/(N-1);

% Initial Value Specification " From where to start "
y0 = 0.5; 

% Identifying variable "t" as a function for time
t(1) = a; % In matlab the zero index is replaced with index 1 

% identifying the variable "y" in matlab 
y = zeros(1,numel(t));

% Implemnet the Modified Euler Methods 
for j = 1:N-1
    t(j+1) = t(j)+h;
    y(j+1)=y(j)+(h)*f(t(j)+h,y(j)+h*f(t(j),y(j)));
    fprintf('%6.2f %12.8f\n',t(j),y(j));
end

plot(t,y,'k-o');
hold on 

% The actual solution 

% The actual solution 
syms tt
syms tt
yy(tt) = piecewise(0<tt<1,5*tt-5*(tt).^2,1<tt<2,-5*(tt).^2);

fplot(yy)

xvalues = linspace(0.1,1.99,N);
yy_values = yy(xvalues);

title(" Backward Euler Method versus Analytical Solutions")
xlabel("Protein Plasma Concentration")
ylabel("Synthesis Rate")
h = legend ('Backward Euler Method','Analytical Solution');
hold on;
grid
% To convert the sym into numerics to evaluate the error
aa = double(yy_values); 
error1=[' error: ' num2str(100*abs((aa(end)-y(end))/aa(end))) '%'];

% Specification of the function 
function f = f(t,y)
if t<1
     f =5-10*t;
else
   f =-10*t;
end
end