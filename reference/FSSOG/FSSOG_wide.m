%For wider Gaussian, low-rank convolution/expansion will be applied on the direction with/without periodicity
function [Pot_wide,realphi,err,t_long,t_Dlong]=FSSOG_wide(M_min,M_max,L,Nx,Ny,R,x,N)

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


%% Step 1: Directly caluculate the charge on Fourier-Chebyshev nodes
t_long=0;
tic;

H_tilde_kspace(1:2*Nx+1,1:2*Ny+1,1:R)=0;
for i=1:2*Nx+1
    for j=1:2*Ny+1
        k_x=kx(i);
        k_y=ky(j);
        kk=sqrt(k_x^2+k_y^2);
        for k=1:R
            z_k=rz(k);
            for t=1:N3
                xt=x(t,1); yt=x(t,2); zt=x(t,3); qt=x(t,4);
                for ell=1:l_range
                    H_tilde_kspace(i,j,k)=H_tilde_kspace(i,j,k)+pi*qt*wl(ell)*(2-4*(z_k-zt)^2/sl(ell)^2+kk^2*sl(ell)^2)*exp(-(z_k-zt)^2/sl(ell)^2)*exp(-sl(ell)^2*kk^2/4)*exp(-1i*(k_x*xt+k_y*yt));
                end
                 %H_tilde_kspace(i,j,k)=H_tilde_kspace(i,j,k)+qt/4*wl*(1/2*Tfunc(z_k,eta)*kk^2-Coefz(z_k,sl,eta))*sl^2*exp(-(z_k-zt)^2/sl^2)*exp(-sl^2*kk^2/4)*exp(-1i*(k_x*xt+k_y*yt));
            end
        end
    end
end   

21



%% Step 2: Convert to Chebyshev series.
B_in(1:2*Nx+1,1:2*Ny+1,1:R)=0;
for k=1:R
  for l=1:R
     tempB(1:2*Nx+1,1:2*Ny+1)= H_tilde_kspace(:,:,l);
     B_in(:,:,k)=B_in(:,:,k)+2/R.*tempB.*ChebPoly(k-1,rz(l),domain);
  end
end

22

t_long=t_long+toc;
%% Step 3: SolEqn
coef_Phi(1:2*Nx+1,1:2*Ny+1,1:R)=0; %For each k-mode, it stores [phi_0'',phi_1'',...,phi_{R-1}'',phi_0,phi_0']

for i=1:2*Nx+1
    for j=1:2*Ny+1
        if (i==Nx+1) && (j==Nx+1)
            continue;
        end
        k_x=kx(i);
        k_y=ky(j);
        kk=sqrt(k_x^2+k_y^2);
        
        
        mu=0; nu=-kk^2;  scale=domain; 
        A(1:R,1:R+2)=0;
        A(1,1)=1;
        A(2,1)=mu/2; A(2,2)=1+nu/8*scale^2; A(2,3)=-mu/2; A(2,4)=-nu/8*scale^2;
        for s=3:R
            k=s-1;
            A(s,s-2)=nu/(2*k*(2*(k-1)))*scale^2;
            A(s,s-1)=mu/(2*k);
            A(s,s)=1-nu/(2*k)*(1/(2*(k-1))+1/(2*(k+1)))*scale^2;
            A(s,s+1)=-mu/(2*k);
            A(s,s+2)=nu/(2*k*(2*(k+1)))*scale^2;
        end

        
        B(1:R,1:2)=0;
        B(1,1)=2*nu; B(1,2)=2*mu; B(2,2)=nu*scale;
        
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
        
        
        alpha=0; beta=0;
        for t=1:N3
            for ell=1:l_range
                alpha=alpha+pi*wl(ell)*sl(ell)^2*exp(-sl(ell)^2*kk^2/4)*x(t,4)*exp(-1i*(k_x*x(t,1)+k_y*x(t,2)))*exp(-(-Lz/2-x(t,3))^2/sl(ell)^2);
                beta=beta+pi*wl(ell)*sl(ell)^2*exp(-sl(ell)^2*kk^2/4)*x(t,4)*exp(-1i*(k_x*x(t,1)+k_y*x(t,2)))*exp(-(Lz/2-x(t,3))^2/sl(ell)^2);
            end
        end

        Divisor(1:R+2,1:R+2)=0;
        Divisor(1:R,1:R)=A(1:R,1:R);
        Divisor(1:R,R+1:R+2)=B(1:R,1:2);
        Divisor(R+1:R+2,1:R)=C(1:2,1:R);
        Divisor(R+1:R+2,R+1:R+2)=D(1:2,1:2);
        
        
        tic
        check=Divisor;
        RHS(1:R+2,1)=0;
        RHS(1:R,1)=-B_in(i,j,:);
        RHS(R+1,1)=alpha;
        RHS(R+2,1)=beta;
        PhiPP(1:R+2)=Divisor\RHS;                                                                                                                                                               
                     
        
        % Calculate the original coefficient
        SD(1:R+2)=0;
        SD(1:R)=PhiPP(1:R);
        coef_Phi(i,j,1)=2*PhiPP(R+1);
        coef_Phi(i,j,2)=1/8*scale^2*(SD(2)-SD(4))+scale*PhiPP(R+2);
        for s=3:R
            k=s-1;
            coef_Phi(i,j,s)=1/(2*k)*scale^2*(1/(2*(k-1))*(SD(s-2)-SD(s))-1/(2*(k+1))*(SD(s)-SD(s+2)));  
        end
        t_long=t_long+toc;
    end
end
23


%% Step 4:Gathering
tic;
phi(1:N3)=0;

for t=1:N3
    xt=x(t,1); yt=x(t,2); zt=x(t,3); qt=x(t,4);
    for i=1:2*Nx+1
        k_x=kx(i);
        for j=1:2*Ny+1
            if ((i==Nx+1) && (j==Ny+1))
                continue;
            end
            k_y=ky(j);
            for k=1:R
                if (k==1)
                    flag=1/2;
                else
                    flag=1;
                end
                phi(t)=phi(t)+1/(Lx*Ly)*flag*coef_Phi(i,j,k)*ChebPoly(k-1,zt,domain)*exp(1i*(k_x*xt+k_y*yt));
            end
        end
    end
end
t_long=t_long+toc;

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
24


t_long=t_long+toc;
%% Verification

tic;

realphi(1:N3)=0;
% cutoff=40;
% 
% for ell=1:l_range
% Kmax=12/sl(ell);
% h=2*Kmax/cutoff;
% Kz=-Kmax:h:Kmax;
% siumsum=0;
% for k1=-16:16
%     k_x=k1*2*pi/Lx;
%     for k2=-16:16
%         k_y=k2*2*pi/Ly;
%         for k3=1:(cutoff+1)
%             k_z=Kz(k3);
%             %if(abs(k1)<=0&&abs(k2)<=0)
%             %    [k_x k_y k_z exp(-1*sl^2*(k_x^2+k_y^2+k_z^2)/4)]
%             %end
%             for i=1:N3
%                 for j=1:N3
%                 
%                        realphi(i)=realphi(i)+(pi^(1/2)/(2*Lx*Ly))*wl(ell)*sl(ell)^3*h*x(j,4)*exp(-1*sl(ell)^2*(k_x^2+k_y^2+k_z^2)/4)*exp(1i*((x(i,1)-x(j,1))*k_x+(x(i,2)-x(j,2))*k_y+(x(i,3)-x(j,3))*k_z));
%                   
%                 end
%             end
%         end
%     end
% end
% end
cutoff=30;

for i=1:N3
    xi=x(i,1); yi=x(i,2); zi=x(i,3); qi=x(i,4);
for j=1:N3
    xj=x(j,1); yj=x(j,2); zj=x(j,3); qj=x(j,4);
    for ell=1:l_range
        for k1=-cutoff:cutoff
            for k2=-cutoff:cutoff
                realphi(i)=realphi(i)+pi/L^2*qj*wl(ell)*sl(ell)^2*exp(-(zi-zj)^2/sl(ell)^2)*exp(-sl(ell)^2*(2*pi/L)^2*(k1^2+k2^2)/4)*exp(1i*(2*pi/L)*(k1*(xi-xj)+k2*(yi-yj)));
            end
        end
    end
end
end


error=(phi-realphi)';
err=max(abs(error));
25

t_Dlong=toc;


