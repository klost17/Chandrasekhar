%% Chandrasekhar limit
%% Matlab codes used to find the numerical solutions
% Justo, Alexandre   Navarro, Pablo   Ruíz, José Javier   Nigorra, Joan
%% FINDING THE LIMIT MASS

clear all
close all

% We solve the system. How z can not be 0, we input that the RK4 method starts integrating in z=10^-6
[zvec,yvec,h]=climit(10^(-5),30,10000,1,0);
%We find for which value of z the solution Y(z) vanishes. Knowing that at
%some point the solution cross the x-axis: 
avec=[find(yvec<0)];
m=avec(1);
if abs(yvec(m-1))<=abs(yvec(m));
    zo=zvec(m-1);
else
    zo=zvec(m);
end
%zo=value for which Y(zo)=0;

zfvec=zvec(1:m); yfvec=yvec(1:m);
gvec=yfvec.^3.*zfvec.^2;


%We calculate the integral as the sum of the areas of the trapezoids. Each
%one with area equal to h/2*(gvec(i)+gvec(i+1)).Notice that we only take
%gvec(1) and gvec(end) once, and the inner components twice. Then:
M=h*(gvec(1)/2 +sum(gvec(2:end-1))+gvec(end)/2);
%The mass limit (Chandrasekhar limit) is:
Chandrasekharlimit=3.678*10^30/(1.989*10^30)*sqrt(3*pi)/8*M 

%RUNGE-KUTTA 4
%{
function [zvec,yvec,h]=climit(zo,zf,n,yo,uo)
Y=[yo;uo];
h=(zf-zo)/n;
zvec=[zo];
yvec=[yo];
for i=1:n
g1=h*funcionderivada(zvec(end),Y);
g2=h*funcionderivada(zvec(end)+h/2,Y+0.5*g1);
g3=h*funcionderivada(zvec(end)+h/2,Y+0.5*g2);
g4=h*funcionderivada(zvec(end)+h,Y+g3);
zvec=[zvec zvec(end)+h];
Y=Y+1/6*g1+1/3*g2+1/3*g3+1/6*g4;
yvec=[yvec Y(1)];
end
plot(zvec,yvec);grid on
figure(1)
xlabel('z');ylabel('Y(z)');hold on;
end

function dY=funcionderivada(z,Y);
y = Y(1);
u = Y(2);

%The system we have to solve is:
dx=u;
du =-y^3-2/z*u; 
dY=[dx;du];
end
%}

%% FINDING THE CHANDRASEKHAR LIMIT

clear all
close all
% Now we are going to study how mass and radius depend on Xc 
Mcvec=[]; Rvec=[];xcvec=[];
for xc=[10^-2:0.01:1 1:1:50]
[zvec,yvec,h]=densities(10^-5,7,10000,1,0,xc);
%As we can see in the system we have to solve, in the du expression appears
%Y(z) in the denominator, so when the solution vanishes we will have a
%problem. Sometimes the solution goes suddenly to -inifint and sometimes to
%+ infinity without crossing the x-axis, consequently we can not use the method we
%used before to find zo. 
t=abs(yvec);
b=min(t);
m=find(t==b);
R=zvec(m)*3854283.3/(xc*6.378*10^6); %in units of Earth radii 
zfvec=zvec(1:m); yfvec=yvec(1:m);
gvec=yfvec.^3.*zfvec.^2;
M=h*(gvec(1)/2 +sum(gvec(2:end-1))+gvec(end)/2);
Mc=sqrt(3*pi)/8*3.678*10^(30)*M/(1.989*10^30); %in solar masses
Mcvec=[Mcvec Mc];
Rvec=[Rvec R];
xcvec=[xcvec xc];
end
mlimit=Mcvec(end);
l=length(xcvec);
mlimitvec=mlimit*linspace(1,1,l);

figure(1)
plot(xcvec,Rvec),grid on
 xlabel('Value of Xc')
 ylabel('White Dwarf Radius (in units of Earth radii)')
figure(2)
plot(xcvec,Mcvec),grid on,hold on
plot(xcvec,mlimitvec,'r--'),grid on
 xlabel('Value of Xc')
 ylabel('White Dwarf Mass (in solar masses)')
figure(3)
plot(Mcvec,Rvec),grid on
 xlabel('White Dwarf Mass (in solar masses)')
 ylabel('White Dwarf Radius (in units of Earth radii)')

% RUNGE-KUTTA 4
%{
function [zvec,yvec,h]=densities(zo,zf,n,yo,uo,xc)
Y=[yo;uo];
h=(zf-zo)/n;
zvec=[zo];
yvec=[yo];
for i=1:n
g1=h*funcionderivada(zvec(end),Y,xc);
g2=h*funcionderivada(zvec(end)+h/2,Y+0.5*g1,xc);
g3=h*funcionderivada(zvec(end)+h/2,Y+0.5*g2,xc);
g4=h*funcionderivada(zvec(end)+h,Y+g3,xc);
zvec=[zvec zvec(end)+h];
Y=Y+1/6*g1+1/3*g2+1/3*g3+1/6*g4;
yvec=[yvec Y(1)];
end
%plot(zvec,yvec);grid on
%figure(1)
%xlabel('z');ylabel('Y(z)');hold on;
end

function dY=funcionderivada(z,Y,xc);
y = Y(1);
u = Y(2);

%The system we have to solve is:
dx=u;
du =-2*u/z-u^2/(y*(1+xc^2*y^2))-y^2*sqrt(1+xc^2*y^2)/xc; %*** Mirar abajo
dY=[dx;du];
end
%}
