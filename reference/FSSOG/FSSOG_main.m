clear

b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
M=25;
L=100;
Lx=L; Ly=L; Lz=L;
Nx=64; Ny=64; Nz=64; 
NLx=16; NLy=16; R=30;
wsearch=Nx/4;
beta=70;
h=L/Nx;
%wsearch=8; %Interpolation num, 5-grids for each direction.
width=wsearch*h; %Interpolation distance


N=5;
N3=N^3;
% x(1,1:4)=[-Lx/3,Ly/3,2,1]; 
% x(2,1:4)=[Lx/3,Ly/3,-2,-2];
% x(3,1:4)=[Lx/3,-Ly/3,2,-3];
% x(4,1:4)=[-Lx/3,Ly/3,-2,-4];
% x(5,1:4)=[Lx/3,-Ly/3,-2,4];
% x(6,1:4)=[-Lx/3,-Ly/3,2,2];
% x(7,1:4)=[-Lx/3,Ly/3,-2,3];
% x(8,1:4)=[-Lx/3,-Ly/3,10,-1];

x_temp=load('x.mat','x');
xx_temp=struct2cell(x_temp);
x=cell2mat(xx_temp);

% %Randomly generate particle information
% x(1:N3,1:4)=0;
% for i=1:N3
%     x(i,1)=Lx*(rand-0.5);
%     x(i,2)=Ly*(rand-0.5);
%     x(i,3)=Lz*(rand-0.5);
%     x(i,4)=2*randn;
% end
% ave_summa=sum(x(:,4))./N3; %Charge Neutrality
% x(:,4)=x(:,4)-ave_summa;

mm=0:1:M;
l_range=length(mm);
wl=(2*log(b))./sqrt(2*pi*sigma^2).*(1./b.^mm);
sl=sqrt(2).*b.^mm*sigma;  


M_critical=0;
ratio=1/5;% The criterion of mid-range and long-range Gaussian;

for ell=1:l_range
    if (sl(ell)<=ratio*L) && (sl(ell+1)>ratio*L)
        M_critical=ell-1;
        break;
    end
end

%% Calculate Potential

% acc_wide_store(1:10)=0;
% for mm=0:9
% mm
% [Pot_wide,realphi_wide,acc_wide,t_long,t_Dlong]=FSSOG_wide(mm,mm,L,NLx,NLy,R,x,N);
% acc_wide_store(mm+1)=acc_wide;
% end
% acc_wide_store
% save('acc_wide.mat','acc_wide_store')

[Pot_mid,realphi_mid,acc_mid,t_mid,t_Dmid]=FSSOG_mid(0,M_critical,L,Nx,Ny,Nz,wsearch,x,N,beta);
[Pot_wide,realphi_wide,acc_wide,t_long,t_Dlong]=FSSOG_wide(M_critical+1,M,L,NLx,NLy,R,x,N);
% [Pot_wide,realphi_wide,acc_wide,t_long,t_Dlong]=FSSOG_wide(9,9,L,NLx,NLy,R,x,N);
Phi=Pot_mid+Pot_wide;

%% Verification
realphi(1:N3)=realphi_mid+realphi_wide;
% cutoff=30;
% 
% for i=1:N3
%     xi=x(i,1); yi=x(i,2); zi=x(i,3); qi=x(i,4);
% for j=1:N3
%     xj=x(j,1); yj=x(j,2); zj=x(j,3); qj=x(j,4);
%     for ell=1:l_range
%         for k1=-cutoff:cutoff
%             for k2=-cutoff:cutoff
%                 realphi(i)=realphi(i)+pi/L^2*qj*wl(ell)*sl(ell)^2*exp(-(zi-zj)^2/sl(ell)^2)*exp(-sl(ell)^2*(2*pi/L)^2*(k1^2+k2^2)/4)*exp(1i*(2*pi/L)*(k1*(xi-xj)+k2*(yi-yj)));
%             end
%         end
%     end
% end
% end


% realphi_2(1:N3)=0;
% scale=0.1;
% cutoff=20;
% for k1=-cutoff:cutoff
%     for k2=-cutoff:cutoff
%         for k3=-cutoff/scale:cutoff/scale
%             k_x=k1*2*pi/L;
%             k_y=k2*2*pi/L;
%             k_z=k3*2*pi/(L+2*width)*scale;
%             for i=1:N3
%                 for j=1:N3
%                     for ell=1:M_critical+1
%                         realphi_2(i)=realphi_2(i)+(pi^(3/2)/(L^2*(L+2*width)))*scale*wl(ell)*sl(ell)^3*x(j,4)*exp(-1*sl(ell)^2*(k_x^2+k_y^2+k_z^2)/4)*exp(1i*((x(i,1)-x(j,1))*k_x+(x(i,2)-x(j,2))*k_y+(x(i,3)-x(j,3))*k_z));          
%                     end
%                 end
%             end
%         end
%     end
% end
% 


err=max(abs(Phi-realphi))

