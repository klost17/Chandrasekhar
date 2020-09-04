%Encontrar masa límite
[zvec,yvec,h]=climit(10^(-6),30,10000,1,0);
avec=[find(yvec<0)];
m=avec(1);
if abs(yvec(m-1))<=abs(yvec(m))
    zo=zvec(m-1)
else
    zo=zvec(m)
end
%zo=valor para el cual y(zo)=0
zfvec=zvec(1:m); yfvec=yvec(1:m);
gvec=yfvec.^3.*zfvec.^2;
figure(2)
plot(zfvec,gvec),grid on; xlabel('z');ylabel('g(z)');

%calculamos integral sumando áreas de trapecios
M=h*(gvec(1)/2 +sum(gvec(2:end-1))+gvec(end)/2);

Mc=sqrt(3*pi)/8*M % *Nh^3*mh