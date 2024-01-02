%For wider Gaussian, low-rank convolution/expansion will be applied on the direction with/without periodicity
function [Pot_wide,err,t_long,t_Dlong]=FSSOG_wide_opt(M_min,M_max,L,Nx,Ny,R,x,N)

%%  Initialization
M=M_min:1:M_max;
l_range=length(M);
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
wl=(2*log(b))./sqrt(2*pi*sigma^2).*(1./b.^M);
sl=sqrt(2).*b.^M*sigma;  


%Particle information(z-direction is free):(x,y,z,q)
Lx=L; Ly=L; Lz=L;

% N=5;
N3=N^3;
% x(1,1:4)=[-Lx/3,Ly/3,2,1]; 
% x(2,1:4)=[Lx/3,Ly/3,-2,-2];
% x(3,1:4)=[Lx/3,-Ly/3,2,-3];
% x(4,1:4)=[-Lx/3,Ly/3,-2,-4];
% x(5,1:4)=[Lx/3,-Ly/3,-2,4];
% x(6,1:4)=[-Lx/3,-Ly/3,2,2];
% x(7,1:4)=[-Lx/3,Ly/3,-2,3];
% x(8,1:4)=[-Lx/3,-Ly/3,10,-1];

%Method setting
%Nx=16; Ny=16; %number of grid points in xy-direction (one-half)
%R=10; %number of grid points in z-directions
hx=Lx/(2*Nx); 
hy=Ly/(2*Ny); 

% rx(1:Nx)=linspace(-Lx/2,Lx/2-hx,Nx);
% ry(1:Ny)=linspace(-Ly/2,Ly/2-hy,Ny); 
kx(1:2*Nx+1)=linspace(-2*pi/Lx*Nx,2*pi/Lx*Nx,2*Nx+1);
ky(1:2*Ny+1)=linspace(-2*pi/Ly*Ny,2*pi/Ly*Ny,2*Ny+1);
rz(1:R)=0;
domain=Lz/2;


for i=1:R
    rz(i)=-domain*cos((2*i-1)/(2*R)*pi); 
end

kk_temp(1:2*Nx+1,1:2*Ny+1)=sqrt(kx'.^2+ky.^2);
k_x_temp(1:2*Nx+1,1:2*Ny+1)=repmat(kx',1,2*Ny+1);
k_y_temp(1:2*Nx+1,1:2*Ny+1)=repmat(ky,2*Nx+1,1);

%% Step 1: Directly caluculate the charge on Fourier-Chebyshev nodes
t_long=0;
H_tilde_kspace(1:2*Nx+1,1:2*Ny+1,1:R)=0;
A2=zeros(2*Nx+1,2*Nx+1,l_range);
for ell=1:l_range
    A2(:,:,ell)=exp(-sl(ell)^2.*kk_temp.^2./4);
end

tic
for t=1:N3
    xt=x(t,1); yt=x(t,2); zt=x(t,3); qt=x(t,4);
    A1=exp(-1i.*(k_x_temp.*xt+k_y_temp.*yt));

    A3(1:2*Nx+1,1:2*Ny+1)=0;
    for ell=1:l_range
        A3=A1.*A2(:,:,ell);
        temp_H(1:2*Nx+1,1:2*Ny+1,1:R)=0;
        for k=1:R
            z_k=rz(k);                        
            temp_H(:,:,k)=temp_H(:,:,k)+pi*qt*wl(ell).*(2-4*(z_k-zt)^2/sl(ell)^2+kk_temp.^2.*sl(ell)^2).*exp(-(z_k-zt)^2/sl(ell)^2);
        end
        H_tilde_kspace=H_tilde_kspace+temp_H.*A3;
    end
end
toc
%% Step 2: Convert to Chebyshev series.
tic
B_in(1:2*Nx+1,1:2*Ny+1,1:R)=0;
for k=1:R
  for l=1:R
     tempB(1:2*Nx+1,1:2*Ny+1)= H_tilde_kspace(:,:,l);
     B_in(:,:,k)=B_in(:,:,k)+2/R.*tempB.*ChebPoly(k-1,rz(l),domain);
  end
end

toc
%t_long=t_long+toc;
%% Step 3: SolEqn
coef_Phi(1:2*Nx+1,1:2*Ny+1,1:R)=0; %For each k-mode, it stores [phi_0'',phi_1'',...,phi_{R-1}'',phi_0,phi_0']

%Basic matrix construction
scale=domain; 
A(1:R,1:R+2)=0;
%A(1,1)=1;
%A(2,2)=1+nu/8*scale^2; A(2,4)=-nu/8*scale^2;

A(2,2)=1/8*scale^2; A(2,4)=-1/8*scale^2;
for s=3:R
    k=s-1;
    A(s,s-2)=1/(2*k*(2*(k-1)))*scale^2;
    A(s,s)=-1/(2*k)*(1/(2*(k-1))+1/(2*(k+1)))*scale^2;
    A(s,s+2)=1/(2*k*(2*(k+1)))*scale^2;
end

        
B(1:R,1:2)=0;
B(1,1)=2; B(2,2)=scale;
        
%A and B with nu=-kk^2, C and D without nu.
C(1:2,1:R+2)=0; 
for s=3:R
    k=s-1;
    C(1,s-2)=C(1,s-2)+scale^2*(1/(2*k))*(1/(2*(k-1)))*(-1)^k;
    C(1,s)=C(1,s)-scale^2*((1/(2*k))*(1/(2*(k-1)))+(1/(2*k))*(1/(2*(k+1))))*(-1)^k;
    C(1,s+2)=C(1,s+2)+scale^2*(1/(2*k))*(1/(2*(k+1)))*(-1)^k;
    
    C(2,s-2)=C(2,s-2)+scale^2*(1/(2*k))*(1/(2*(k-1)));
    C(2,s)=C(2,s)-scale^2*((1/(2*k))*(1/(2*(k-1)))+(1/(2*k))*(1/(2*(k+1))));
    C(2,s+2)=C(2,s+2)+scale^2*(1/(2*k))*(1/(2*(k+1)));
end
C(1,2)=C(1,2)+1/8*scale^2*(-1);
C(1,4)=C(1,4)-1/8*scale^2*(-1);
        
C(2,2)=C(2,2)+1/8*scale^2;
C(2,4)=C(2,4)-1/8*scale^2;
        
        
D(1:2,1:2)=0;
D(1,1)=1;
D(1,2)=-scale;
D(2,1)=1;
D(2,2)=scale;

tic
alpha(1:2*Nx+1,1:2*Ny+1)=0;
beta(1:2*Nx+1,1:2*Ny+1)=0;
for t=1:N3
    for ell=1:l_range
        alpha=alpha+pi*wl(ell)*sl(ell)^2.*exp(-sl(ell)^2.*kk_temp.^2./4).*x(t,4).*exp(-1i.*(k_x_temp.*x(t,1)+k_y_temp.*x(t,2))).*exp(-(-Lz/2-x(t,3))^2/sl(ell)^2);
        beta=beta+pi*wl(ell)*sl(ell)^2.*exp(-sl(ell)^2.*kk_temp.^2./4).*x(t,4).*exp(-1i.*(k_x_temp.*x(t,1)+k_y_temp.*x(t,2))).*exp(-(Lz/2-x(t,3))^2/sl(ell)^2);
    end
end
toc
% t_long=t_long+toc;


%Store the inverse
Divisor_inv(1:R+2,1:R+2,1:2*Nx+1,1:2*Ny+1)=0;
for i=1:2*Nx+1
    for j=1:2*Ny+1
        if (i==Nx+1) && (j==Nx+1)
            continue;
        end
        k_x=kx(i);
        k_y=ky(j);
        kk=sqrt(k_x^2+k_y^2);       
        
        nu=-kk^2;              

        Divisor(1:R+2,1:R+2)=0;
        Divisor(1:R,1:R)=nu.*A(1:R,1:R)+eye(R);
        Divisor(1:R,R+1:R+2)=nu.*B(1:R,1:2);
        Divisor(R+1:R+2,1:R)=C(1:2,1:R);
        Divisor(R+1:R+2,R+1:R+2)=D(1:2,1:2);    
        
        Divisor_inv(:,:,i,j)=inv(Divisor);
    end
end

%SolEqn
tic
for i=1:2*Nx+1
    for j=1:2*Ny+1
        if (i==Nx+1) && (j==Nx+1)
            continue;
        end
        RHS(1:R+2,1)=0;
        RHS(1:R,1)=-B_in(i,j,:);
        RHS(R+1,1)=alpha(i,j);
        RHS(R+2,1)=beta(i,j);
        
        D_inv(1:R+2,1:R+2)=Divisor_inv(:,:,i,j);
        PhiPP(1:R+2)=D_inv*RHS;                                                                                                                                                                                    
        
        % Calculate the original coefficient
        SD(1:R+2)=0;
        SD(1:R)=PhiPP(1:R);
        coef_Phi(i,j,1)=2*PhiPP(R+1);
        coef_Phi(i,j,2)=1/8*scale^2*(SD(2)-SD(4))+scale*PhiPP(R+2);
        for s=3:R
            k=s-1;
            coef_Phi(i,j,s)=1/(2*k)*scale^2*(1/(2*(k-1))*(SD(s-2)-SD(s))-1/(2*(k+1))*(SD(s)-SD(s+2)));  
        end
    end
end
toc
% t_long=t_long+toc;


%% Step 4:Gathering
tic
phi(1:N3)=0;

for t=1:N3
    xt=x(t,1); yt=x(t,2); zt=x(t,3); 
    for k=1:R
        if (k==1)
           flag=1/2;
        else
           flag=1;
        end
        phi(t)=phi(t)+sum(sum(1/(Lx*Ly)*flag.*coef_Phi(:,:,k).*ChebPoly(k-1,zt,domain).*exp(1i.*(k_x_temp.*xt+k_y_temp.*yt))));
    end
end


toc
%t_long=t_long+toc;

%Zero-mode FCT
Coef=0;
for ell=1:l_range
    Coef=Coef+pi/(Lx*Ly)*wl(ell)*sl(ell)^2.*coef_precompute(sl(ell)^2/(2*domain)^2,2*domain,R);
end
source(1:N3)=x(:,3);
charge(1:N3)=x(:,4);


tic;
mu(1:R)=0;
for l=1:R
  for j=1:N3
      mu(l)=mu(l)+ChebPoly(l-1,source(j),domain)*charge(j);
  end
end

lambda(1:R)=0;
for k=1:R
    for l=1:R+1-k
        lambda(k)=lambda(k)+Coef(k,l)*mu(l);
    end
end

temp_phi(1:N3)=0;
for i=1:N3
    for k=1:R
        temp_phi(i)=temp_phi(i)+lambda(k)*ChebPoly(k-1,source(i),domain);
    end
end
    
phi=phi+temp_phi;

Pot_wide=phi;
toc

%t_long=t_long+toc;
%% Verification

[real_phi,t_Dlong]=Wide_Direct(M_min,M_max,L,x,N);

error=(phi-real_phi)';
err=max(abs(error));



