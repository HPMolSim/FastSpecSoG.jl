MM=[3 4 10 20];
L=100;
Nx=32; Ny=32; Nz=32;
wsearch=Nx/4;
N3=125;
x_temp=load('x.mat','x');
xx_temp=struct2cell(x_temp);
x=cell2mat(xx_temp);

pad_rat(1:100)=0;
for i=1:100
    pad_rat(i)=i;
end


Acc(1:4,1:100)=0;

for t=1:4
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
wl=(2*log(b))/sqrt(2*pi*sigma^2)*(1/b^MM(t));
sl=sqrt(2)*b^MM(t)*sigma;
    
realphi(1:N3)=0;
cutoff=30;
t

for i=1:N3
    xi=x(i,1); yi=x(i,2); zi=x(i,3); qi=x(i,4);
for j=1:N3
    xj=x(j,1); yj=x(j,2); zj=x(j,3); qj=x(j,4);
        for k1=-cutoff:cutoff
            for k2=-cutoff:cutoff
                realphi(i)=realphi(i)+pi/L^2*qj*wl*sl^2*exp(-(zi-zj)^2/sl^2)*exp(-sl^2*(2*pi/L)^2*(k1^2+k2^2)/4)*exp(1i*(2*pi/L)*(k1*(xi-xj)+k2*(yi-yj)));
            end
        end
end
end
    
    for i=1:100
        i
        Acc(t,i)=FSSOG_mid_padding(MM(t),L,Nx,Ny,Nz,wsearch,pad_rat(i),realphi);
    end
end
save('Acc.mat','Acc');

plot(pad_rat,Acc(1,:),'r');
hold on;
plot(pad_rat,Acc(2,:),'b');
hold on;
plot(pad_rat,Acc(3,:),'g');
hold on;
plot(pad_rat,Acc(4,:),'p');
hold on;
legend("M=3","M=4","M=10","M=20");
