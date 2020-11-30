function [  ] = Two_Dimensinal_Differential_Equation_Of_Second_Order1( )

clc
clear all
tstart = 0;   % initial time
tend   = 1000; % final time
N = 100000;  % # of sampling points of solution
energy=zeros(N,1);

tspan = linspace(tstart,tend,N);
%
% define initial conditions and error tolerance
x0=0.1;
xdot0=0.3;
y0=0.1;
ydot0=-0.021;
z0=0;
zdot0=-1;

global g k  m 
m=1;
k=0.001;
g=9.81;


options = odeset('RelTol',1.0e-13,'AbsTol',1.0e-14);
%
% call ode45
[T,outdata] = ode45(@DifferentialFunction,tspan,[x0,xdot0,y0,ydot0,z0,zdot0],options);
x=outdata(:,1);
xdot=outdata(:,2);
y=outdata(:,3);
ydot=outdata(:,4);
z=outdata(:,5);
zdot=outdata(:,6);

 r=sqrt(x.^2+y.^2);
%   z=-k./r;%instead of z from ode45
%   energy =0.5*(xdot.^2+ydot.^2+zdot.^2)-g*z;
% plot(energy);


 hold on
plot3(x,y,z);
 xlabel('x');
 ylabel('y');
 zlabel('z');
 grid on;
 
%  plx(1)=x(1);plx(2)=x(1);ply(1)=y(1);ply(2)=y(1);plz(1)=z(1);plz(2)=z(1);
%  h=plot3(plx,ply,plz,'-co','MarkerEdgeColor','k','MarkerFaceColor','b','DisplayName','satellite','MarkerSize',10);
%  set(h,'XData',plx,'YData',ply,'ZData',plz);
%  for k=2:20:N
%      plx(1)=x(k);plx(2)=x(k);ply(1)=y(k);ply(2)=y(k);plz(1)=z(k);plz(2)=z(k);
%      set(h,'XData',plx,'YData',ply,'ZData',plz);
%      drawnow
%  end 
 hold off

end
%% main function
 
function outputarray = DifferentialFunction(t,indata);
%

global g k  m

c=k;
x=indata(1);xdot=indata(2);y=indata(3);ydot=indata(4); zdot=indata(6);
r=sqrt(x^2+y^2);

help5=y^2-2*x^2;
help6=x^2-2*y^2;
help75=6*x*y*xdot*ydot;
j1=xdot^2*help5+ydot^2*help6-help75;
m1=c/r^5;
s1=x^3+x*y^2;
a1=y^3+x^2*y;
b1=c*m*y/r^3;
n1=c*m*x/r^3;
% xdotdot=-(g-m(xdotdot*s+ydotdot*a+j)*n
% ydotdot=-(g-m(xdotdot*s+ydotdot*a+j)*B
xdotdot=(j1*m1*n1-g*n1)/(-b1*m1*a1-m1*n1*s1+1);
ydotdot=(j1*m1*b1-g*b1)/(-b1*m1*a1-m1*n1*s1+1);
zdotdot=-m1*(xdotdot*s1+ydotdot*a1*j1);
outputarray= [xdot; xdotdot;ydot;ydotdot;zdot;zdotdot]; 
end