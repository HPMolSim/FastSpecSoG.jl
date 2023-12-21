
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
M=25;
L=100;

xx=linspace(-1000,1000,2001);
M_test=[4,8,12];
y1(1:2001)=0;
y2(1:2001)=0;
y3(1:2001)=0;

wl=(2*log(b))./sqrt(2*pi*sigma^2).*(1./b.^M_test);
sl=sqrt(2).*b.^M_test*sigma;  

for i=1:2001
    y1(i)=exp(-xx(i)^2/sl(1)^2);
    y2(i)=exp(-xx(i)^2/sl(2)^2);
    y3(i)=exp(-xx(i)^2/sl(3)^2);
end

plot(xx,y1,'r');
hold on;
plot(xx,y2,'b');
hold on;
plot(xx,y3,'m');
hold on;
plot([L,L],[0,1],'g--');
hold on;
plot([-L,-L],[0,1],'g--');
hold on;

legend("M=4","M=8","M=12","Box length");