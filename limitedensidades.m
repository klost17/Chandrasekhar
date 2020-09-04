
Mcvec=[]; Rvec=[];xcvec=[];
for xc=[10^-2:0.01:1 1:1:400]
[zvec,yvec,h]=densities(10^-5,7,10000,1,0,xc);
t=abs(yvec);
b=min(t);
m=find(t==b);
R=zvec(m)*3854283.3/(xc*6.378*10^6); %en radios terrestres 
zfvec=zvec(1:m); yfvec=yvec(1:m);
gvec=yfvec.^3.*zfvec.^2;
M=h*(gvec(1)/2 +sum(gvec(2:end-1))+gvec(end)/2);
Mc=sqrt(3*pi)/8*3.678*10^(30)*M/(1.989*10^30); %en masas solares
Mcvec=[Mcvec Mc];
Rvec=[Rvec R];
xcvec=[xcvec xc];
end
plot(Mcvec,Rvec),grid on