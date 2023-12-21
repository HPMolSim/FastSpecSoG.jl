function err=Test(Bcut,w0,base,sigma)
L=20;
N=8;
x(1,1:4)=[L/2,L/2,L/2,1];
x(2,1:4)=[L/2,L/2,-L/2,-1];
x(3,1:4)=[L/2,-L/2,L/2,-1];
x(4,1:4)=[-L/2,L/2,L/2,-1];
x(5,1:4)=[L/2,-L/2,-L/2,1];
x(6,1:4)=[-L/2,-L/2,L/2,1];
x(7,1:4)=[-L/2,L/2,-L/2,1];
x(8,1:4)=[-L/2,-L/2,-L/2,-1];

% x_temp=load('x.mat','x');
% xx_temp=struct2cell(x_temp);
% x=cell2mat(xx_temp);

% x(1,1:4)=[0,0,0,1];
% x(2,1:4)=[1,1,1,-1];
% x(3,1:4)=[-1.5,1,-1.5,1];
% x(4,1:4)=[5,1,1.5,-1];


% Bcut=8;
% w0=0.994446492762232252;
% base=2;
% sigma=5.027010924194599;

source=x;
phi=zeros(N,1);
phi_RBSOG=zeros(N,1);
delx=zeros(N,1);
dely=zeros(N,1);
delz=zeros(N,1);
r=zeros(N,1);
Nmax=30;

% for i=-Nmax:Nmax
%     i
%     source(:,1)=x(:,1)+i*2*L;
%     for j=-Nmax:Nmax
%         source(:,2)=x(:,2)+j*2*L;
%             for m=1:N
%                 delx=x(:,1)-source(m,1);
%                 dely=x(:,2)-source(m,2);
%                 delz=x(:,3)-source(m,3);
%                 r=sqrt(delx.*delx+dely.*dely+delz.*delz);
%                 
%                 for p=1:N
%                     if (r(p)>0.01)
%                         phi(p,1)=phi(p,1)+0.5*x(p,4)*x(m,4).*(-delz(p,1)/r(p)^3);
%                     end
%                 end
%             end
%             %end
%     end
% end
% potentiaL_force=phi(1,1)
% pause(2);

for i=-Nmax:Nmax
    i
    source(:,1)=x(:,1)+i*2*L;
    for j=-Nmax:Nmax
        source(:,2)=x(:,2)+j*2*L; 
            %if (i^2+j^2+k^2<400000)
            for m=1:N
                delx=x(:,1)-source(m,1);
                dely=x(:,2)-source(m,2);
                delz=x(:,3)-source(m,3);
                r=sqrt(delx.*delx+dely.*dely+delz.*delz);
                for p=1:N
                    if(r(p)>=10)
                        phi_RBSOG(p,1)=phi_RBSOG(p,1)+0.5*x(p,4)*x(m,4).*(-delz(p,1)./r(p)^3-SOG(Bcut,w0,base,sigma,r(p),delz(p,1)));%1./r(p);%
                    end
                end
            end
    end
end


err=abs((phi_RBSOG(1,1))/(3.726176847183476e-04));
%base^(2*Bcut-2)*exp(-2*pi*(base^(2*Bcut)-base^2)*(sigma/(2*L))^2)
%% b=2 sigma=5.027010924194599 w0=0.994446492762232252 %%
%[1 2 4 8 16 32 64]
%[0.726271292495564 0.370923663765746 0.087877902122996 -6.676143250049517e-04 -0.006547591468934 -0.00655 -0.006570650201871]

%% b=1.32070036405934420 sigma=2.277149356440992 w0=1.00188914114811980%%
%[1 2 4 8 16 32 64]
%[0.999999999254556 0.999992813074544 0.978631111595203 0.449451135148611 0.048723319898033 5.685446030675484e-04 -5.973827940060706e-08]