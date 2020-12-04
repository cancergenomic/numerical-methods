% Midpoint Method 
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
% t = 0 to t = 1
tspan = [0.01 1];
y0 = 0;
options = odeset('RelTol',1e-9,'AbsTol',1e-8);
[T,Y] = ode45(@f,tspan,y0,options);

% t = 1 to t = 2
tspan0 = [1 2];
y0=0;
[t,y] = ode45(@f,tspan0,y0,options);
t = [T;t(2:end)];  % Remove first value as it will be same as last of previous run
y = [Y;y(2:end,:)];

plot(t,y,'-o') 
hold on;

% The actual solution 

syms tt
yy(tt) = piecewise(0<tt<1,5*tt-5*(tt).^2,1<tt<2,-5*(tt).^2);

fplot(yy)

xvalues = linspace(0.1,1.99,N);
yy_values = yy(xvalues);

title(" Adapted Method versus Analytical Solutions")
xlabel("Protein Plasma Concentration")
ylabel("Synthesis Rate")
h = legend ('Adapted Method','Analytical Solution');
hold on;
grid
% To convert the sym into numerics to evaluate the error
aa = double(yy_values); 
error1=[' error: ' num2str(100*abs((aa(end)-y(end))/aa(end))) '%'];

% Specification of the function 



function f = f(t,y)
if  t<1
    f =-10*(+t);
    
else
     f =5-10*t;  
    
    
end
end