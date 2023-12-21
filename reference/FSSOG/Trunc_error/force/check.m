b1=2;
s1=5.027010924194599;
w1=0.994446492762232252;
w_m1=(pi/2)^(-1)*b1*s1^(-1)*log(b1);

b2=1.62976708826776469;
s2=3.633717409009413;
w2=1.00780697934380681;
w_m2=(pi/2)^(-1)*b2*s2^(-1)*log(b2);

b3=1.48783512395703226;
s3=2.662784519725113;
w3=0.991911705759818;
w_m3=(pi/2)^(-1)*b3*s3^(-1)*log(b3);

b4=1.32070036405934420;
s4=2.277149356440992 ;
w4=1.00188914114811980 ;
w_m4=(pi/2)^(-1)*b4*s4^(-1)*log(b4);

b5=1.21812525709410644;
s5=1.774456369233284;
w5=1.00090146156033341;
w_m5=(pi/2)^(-1)*b5*s5^(-1)*log(b5);

Term=[1 2 3 5 7 10 20 30 50 70 100 200 300];
% F1(1:13,1:1)=0;
% F2(1:13,1:1)=0;
% F3(1:13,1:1)=0;
% F4(1:13,1:1)=0;
% F5(1:13,1:1)=0;
lins=linspace(1,300,300);

limit1(1:300,1:1)=0;
limit2(1:300,1:1)=0;
limit3(1:300,1:1)=0;
limit4(1:300,1:1)=0;
limit5(1:300,1:1)=0;
Q=8;
rc=10;


spconst=0.6;
for i=1:300
%     F1(i,1)=Test(Term(i),w1,b1,s1);
%     F2(i,1)=Test(Term(i),w2,b2,s2);
%     F3(i,1)=Test(Term(i),w3,b3,s3);
%     F4(i,1)=Test(Term(i),w4,b4,s4);
%     F5(i,1)=Test(Term(i),w5,b5,s5);
    limit1(i,1)=(log(b1))^(-3/2)*sqrt(2)*25/3*exp(-pi^2/(2*log(b1)))+2/sqrt(pi)*Q/(b1-1)*log(b1)/(sqrt(2)*s1)*exp(-spconst*b1^lins(i))+w_m1*(sqrt(2)*b1^(-1)*s1)^2*exp(-rc^2/(sqrt(2)*b1^(-1)*s1)^2);
    limit2(i,1)=(log(b2))^(-3/2)*sqrt(2)*25/3*exp(-pi^2/(2*log(b2)))+2/sqrt(pi)*Q/(b2-1)*log(b2)/(sqrt(2)*s2)*exp(-spconst*b2^lins(i))+w_m2*(sqrt(2)*b2^(-1)*s2)^2*exp(-rc^2/(sqrt(2)*b2^(-1)*s2)^2);
    limit3(i,1)=(log(b3))^(-3/2)*sqrt(2)*25/3*exp(-pi^2/(2*log(b3)))+2/sqrt(pi)*Q/(b3-1)*log(b3)/(sqrt(2)*s3)*exp(-spconst*b3^lins(i))+w_m3*(sqrt(2)*b3^(-1)*s3)^2*exp(-rc^2/(sqrt(2)*b3^(-1)*s3)^2);
    limit4(i,1)=(log(b4))^(-3/2)*sqrt(2)*25/3*exp(-pi^2/(2*log(b4)))+2/sqrt(pi)*Q/(b4-1)*log(b4)/(sqrt(2)*s4)*exp(-spconst*b4^lins(i))+w_m4*(sqrt(2)*b4^(-1)*s4)^2*exp(-rc^2/(sqrt(2)*b4^(-1)*s4)^2);
    limit5(i,1)=(log(b5))^(-3/2)*sqrt(2)*25/3*exp(-pi^2/(2*log(b5)))+2/sqrt(pi)*Q/(b5-1)*log(b5)/(sqrt(2)*s5)*exp(-spconst*b5^lins(i))+w_m5*(sqrt(2)*b5^(-1)*s5)^2*exp(-rc^2/(sqrt(2)*b5^(-1)*s5)^2);
    

end

% limit(1:5)=0;
% limit(1)=2*sqrt(2)*exp(-pi^2/(2*log(b1)));
% limit(2)=2*sqrt(2)*exp(-pi^2/(2*log(b2)));
% limit(3)=2*sqrt(2)*exp(-pi^2/(2*log(b3)));
% limit(4)=2*sqrt(2)*exp(-pi^2/(2*log(b4)));
% limit(5)=2*sqrt(2)*exp(-pi^2/(2*log(b5)));

