% The Runge Kutta Method
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

% Implementation of the Runge Kutta 
for k=1:N+1
 t(k+1) = t(k)+h;
 k1 = h*f(t(k),y(k));
 k2 = h*f(t(k)+0.5*h,y(k)+0.5*h*k1);
 k3 = h*f((t(k)+0.5*h),(y(k)+0.5*h*k2));
 k4 = h*f((t(k)+h),(y(k)+h*k3));
 y(k+1) = y(k) + (h/6)*(k1+2*k2+2*k3+k4); 
 fprintf('%6.2f %12.8f\n',t(k),y(k));

end

for i = 4:N
    t(i+1) = t(i)+h;
    y(i+1) = y(i) + (h/24)*(55*f(t(i),y(i))-59*f(t(i-1),y(i-1))+37*f(t(i-2),y(i-2))-9*f(t(i-3),y(i-3)));
    fprintf('%6.2f %12.8f\n',t(i),y(i));
end
    
for j = 1:N
    t(j)=t(j+1);
    y(i)=y(i+1);
end

plot(t,y)

hold on 

% The actual solution 
syms tt
syms tt
yy(tt) = piecewise(0<tt<1,5*tt-5*(tt).^2,1<tt<2,-5*(tt).^2);

fplot(yy)

xvalues = linspace(0.1,1.99,N);
yy_values = yy(xvalues);

title(" Runge Kutta Method versus Analytical Solutions")
xlabel("Protein Plasma Concentration")
ylabel("Synthesis Rate")
h = legend ('Runge Kutta Method','Analytical Solution');
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
