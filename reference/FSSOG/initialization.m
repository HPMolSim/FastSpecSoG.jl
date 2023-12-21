%%  Initialization

N=5;
N3=N^3;
L=100;
x(1:N3,1:4)=0;
for i=1:N3
    x(i,1)=L*(rand-0.5);
    x(i,2)=L*(rand-0.5);
    x(i,3)=L*(rand-0.5);
    x(i,4)=2*randn;
end
ave_summa=sum(x(:,4))./N3; %Charge Neutrality
x(:,4)=x(:,4)-ave_summa;
save('x.mat','x');