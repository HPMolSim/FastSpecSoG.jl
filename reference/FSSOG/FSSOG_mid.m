function  [Pot_mid,realphi,acc,t_mid,t_Dmid]=FSSOG_mid(M_min,M_max,L,Nx,Ny,Nz,wsearch,x,N,beta)



%%  Initialization

%u-series setting——Single Gaussian term
%M=0; %Test Gaussian term
M=M_min:1:M_max;
l_range=length(M);
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
wl=(2*log(b))./sqrt(2*pi*sigma^2).*(1./b.^M);
sl=sqrt(2).*b.^M*sigma;  


%Particle information(z-direction is free):(x,y,z,q)
%L=40;
Lx=L; Ly=L; Lz=L;

% N=4;
N3=N^3;


%Method setting
%Nx=32; Ny=32; Nz=32; %number of grid points in xyz-directions
h=L/Nx;
%wsearch=8; %Interpolation num, 5-grids for each direction.
width=wsearch*h; %Interpolation distance
width_window=width+0.5*h;
%beta=20; %5 10 15 20
rx(1:Nx+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Nx+2*wsearch+2); %With virtual images
ry(1:Ny+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Ny+2*wsearch+2); %With virtual images
rz(1:Nz+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Nz+2*wsearch+2);
% kx(1:Nx)=linspace(-(Nx-1)*pi/L,(Nx-1)*pi/L,Nx);
% ky(1:Ny)=linspace(-(Ny-1)*pi/L,(Ny-1)*pi/L,Ny);
% kz(1:Nz+2*wsearch)=linspace(-(Nz+2*wsearch-1)*pi/L,(Nz+2*wsearch-1)*pi/L,Nz+2*wsearch);
% rx_true(1:Nx)=linspace(-L/2,L/2,Nx);
% ry_true(1:Ny)=linspace(-L/2,L/2,Ny);





%% Step1: GridH. Interpolate each q_i to corresponding girds.
tic

H_temp(1:Nx+2*wsearch+2,1:Ny+2*wsearch+2,1:Nz+2*wsearch+2)=0; 

%Interpolate each particle to corresponding grids
for j=1:N3 
    xj=x(j,1); yj=x(j,2); zj=x(j,3);
    
    % Find the neareast grid num
    xnear=floor((xj+L/2+h+width+0.5*h)/h)+1;
    ynear=floor((yj+L/2+h+width+0.5*h)/h)+1;
    znear=floor((zj+L/2+h+width+0.5*h)/h)+1;
    
    %Interpolate charge
    for x_ind=-1*wsearch:wsearch
        for y_ind=-1*wsearch:wsearch
            for z_ind=-1*wsearch:wsearch
                
                z_interp=z_ind+znear;
                z_dist=abs(rz(z_interp)-zj);
                
                x_interp=x_ind+xnear; 
                x_dist=abs(rx(x_interp)-xj);
                
                y_interp=y_ind+ynear; 
                y_dist=abs(ry(y_interp)-yj);
                
%                 x_dist
%                 y_dist
%                 z_dist
%                 pause(5)
                %Charge Interpolation
                H_temp(x_interp,y_interp,z_interp)= H_temp(x_interp,y_interp,z_interp)+x(j,4)*Wkb(x_dist,width_window,beta)*Wkb(y_dist,width_window,beta)*Wkb(z_dist,width_window,beta);
            end
        end
    end 
end

HH=sum(sum(sum(H_temp)))


%Periodic direction gathering
for i=wsearch+2:Nx+wsearch+1
    for j=wsearch+2:Ny+wsearch+1
        for k=1:Nz+2*wsearch+2
            
            % Four edge
            if (i+Nx<=Nx+2*wsearch+2)
                H_temp(i,j,k)=H_temp(i,j,k)+H_temp(i+Nx,j,k);
            end
            
            if (i-Nx>0)
                H_temp(i,j,k)=H_temp(i,j,k)+H_temp(i-Nx,j,k);
            end
            
            if (j+Ny<=Ny+2*wsearch+2)
                H_temp(i,j,k)=H_temp(i,j,k)+H_temp(i,j+Ny,k);
            end
            
            if (j-Ny>0)
                H_temp(i,j,k)=H_temp(i,j,k)+H_temp(i,j-Ny,k);
            end
            
            %Four corner
            if ((i+Nx<=Nx+2*wsearch+2) && (j+Ny<=Ny+2*wsearch+2)) %Rright-Down
                H_temp(i,j,k)=H_temp(i,j,k)+H_temp(i+Nx,j+Ny,k);
            end
            
            if ((i-Nx>0) && (j+Ny<=Ny+2*wsearch+2)) %Left-Down
                H_temp(i,j,k)=H_temp(i,j,k)+H_temp(i-Nx,j+Ny,k);
            end
            
            if ((i+Nx<=Nx+2*wsearch+2) && (j-Ny>0)) %Rright-Up
                H_temp(i,j,k)=H_temp(i,j,k)+H_temp(i+Nx,j-Ny,k);
            end
            
            if ((i-Nx>0) && (j-Ny>0)) %Left-Up
                H_temp(i,j,k)=H_temp(i,j,k)+H_temp(i-Nx,j-Ny,k);
            end
        end
    end
end


H(1:Nx,1:Ny,1:Nz+2*wsearch+2)=H_temp(wsearch+2:Nx+wsearch+1,wsearch+2:Ny+wsearch+1,1:Nz+2*wsearch+2);

11

%% Step2: HF_RtK. Use FFT to get H_hat
H_kspace(1:Nx,1:Ny,1:Nz+2*wsearch+2)=fftn(H);

12
%% Step3: ScalingH. Take the scaling of eack k-mode
%TdK(1:Nx+2*wsearch,1:Ny+2*wsearch,1:Nz+2*wsearch)=0;
H_tilde_kspace(1:Nx,1:Ny,1:Nz+2*wsearch+2)=0;

% Calculating scaling factor
for i=1:Nx
    for j=1:Ny
        for k=1:Nz+2*wsearch+2
             fx=i-1; fy=j-1; fz=k-1;
            if (fx>ceil(Nx/2))
                fx=fx-Nx;
            end
            if (fy>ceil(Ny/2))
                fy=fy-Ny;
            end
            if (fz>(ceil(Nz/2)+wsearch+1))
                fz=fz-Nz-2*wsearch-2;
            end    
            k_x=fx*2*pi/L;
            k_y=fy*2*pi/L;
            k_z=fz*2*pi/(L+2*width+2*h);
            TdK=0;
            for ell=1:l_range
                TdK=TdK+sqrt(pi)/4*wl(ell)*sl(ell)^3*exp(-1*sl(ell)^2*(k_x^2+k_y^2+k_z^2)/4);
            end
            H_tilde_kspace(i,j,k)=TdK*(FWkb(k_x,width_window,beta)*FWkb(k_y,width_window,beta)*FWkb(k_z,width_window,beta))^(-2)*H_kspace(i,j,k);
        end
    end
end

13
%% Step4:HF_KtR. Use IFFT to get H_tilde
H_tilde(1:Nx,1:Ny,1:Nz+2*wsearch+2)=ifftn(H_tilde_kspace);

14

%% Step5:Numerical Integral for potential
phi(1:N3)=0; %Potential
for j=1:N3 %Contribute the j-th charge
    xj=x(j,1); yj=x(j,2); zj=x(j,3);
    
    % Find the neareast grid num
    xnear=floor((xj+L/2+0.5*h)/h)+1;
    ynear=floor((yj+L/2+0.5*h)/h)+1;
    znear=floor((zj+L/2+h+width+0.5*h)/h)+1;
    
     for x_ind=-1*wsearch:wsearch
        for y_ind=-1*wsearch:wsearch
            for z_ind=-1*wsearch:wsearch
                
                z_interp=z_ind+znear;
                z_dist=abs(rz(z_interp)-zj);
                
                x_temp=x_ind+xnear;
                x_interp=mod(x_temp,Nx);
                if (x_interp==0)
                    x_interp=Nx;
                end
                x_dist=abs((-L/2+(x_temp-1)*h)-xj);
                
                y_temp=y_ind+ynear;
                y_interp=mod(y_temp,Ny);
                if (y_interp==0)
                    y_interp=Ny;
                end
                y_dist=abs((-L/2+(y_temp-1)*h)-yj);
%                 
%                x_dist
%                y_dist
%                z_dist
%                pause(5)
                phi(j)=phi(j)+4*pi*h^3*H_tilde(x_interp,y_interp,z_interp)*Wkb(x_dist,width_window,beta)*Wkb(y_dist,width_window,beta)*Wkb(z_dist,width_window,beta);
            end
        end
     end
end

Pot_mid=phi;

t_mid=toc;

15

%% Verification
tic;

%phi
realphi(1:N3)=0;
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


%realphi
%error=abs(phi-realphi)./abs(realphi);
error=abs(phi-realphi);
acc=max(error);

16

t_Dmid=toc;