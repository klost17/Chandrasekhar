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

dx=u;
du =-y^3-2/z*u; 
dY=[dx;du];
end