function Coef=PKB_coef(h,P,width,beta,nu)
%h:grid size. P:window number; 
%width=width_window=0.5h+Ph+0.5h
%beta=2.5P
%nu=min{9,P/2+2}+1

width_window=width+0.5*h;

C(1:P+2,1:nu)=0;
x(1:nu)=0;
for k=1:nu
    x(k)=-cos(pi*(2*k-1)/(2*nu));
end

Vdm(1:nu,1:nu)=0;
for s=1:nu
   for t=1:nu
       Vdm(s,t)=x(s)^(t-1);
   end
end


%Middle coefficients
for i=2:P+1
   center=-width+(i-3/2)*h;
   scale=1/2*h;
   b(1:nu,1:1)=0;
   for s=1:nu
       b(s,1)=Wkb(x(s)*scale+center,width_window,beta);
   end
   temp=Vdm\b;
   C(i,:)=temp';
end

%i=1
center=-width-1/4*h;
scale=1/4*h;
b(1:nu,1:1)=0;
for s=1:nu
    b(s,1)=Wkb(x(s)*scale+center,width_window,beta);
end
temp=Vdm\b;
C(1,:)=temp';

%i=P+2
center=width+1/4*h;
scale=1/4*h;
b(1:nu,1:1)=0;
for s=1:nu
    b(s,1)=Wkb(x(s)*scale+center,width_window,beta);
end
temp=Vdm\b;
C(P+2,:)=temp';

Coef=C;