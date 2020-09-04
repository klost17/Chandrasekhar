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

dx=u;
du =-2*u/z-u^2/(y*(1+xc^2*y^2))-y^2*sqrt(1+xc^2*y^2)/xc; %*** Mirar abajo
dY=[dx;du];
end

%*** Peta para y==0, hay que vigilar. No vale el método de encontrar el
%primer punto negativo porque a veces peta hacia +infinito sin cruzar el
%eje x y a veces hacia -infinito