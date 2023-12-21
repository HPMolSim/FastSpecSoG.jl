% function  acc=FSSOG_mid_padding(M,L,Nx,Ny,Nz,wsearch,pad_ratio,realphi)
%%  Initialization

% acc(1:10)=0;
% for mm=0:9
%u-series setting——Single Gaussian term
M=13; %Test Gaussian term
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
wl=(2*log(b))/sqrt(2*pi*sigma^2)*(1/b^M);
sl=sqrt(2)*b^M*sigma;



%Particle information(z-direction is free):(x,y,z,q)
L=100;
N=5;
N3=N^3;

x_temp=load('x.mat','x');
xx_temp=struct2cell(x_temp);
x=cell2mat(xx_temp);
% x(1:N3,1:4)=0;
% for i=1:N3
%     x(i,1)=L*(rand-0.5);
%     x(i,2)=L*(rand-0.5);
%     x(i,3)=L*(rand-0.5);
%     x(i,4)=2*randn;
% end
% ave_summa=sum(x(:,4))./N3; %Charge Neutrality
% x(:,4)=x(:,4)-ave_summa;


%Method setting
Nx=16; Ny=16; Nz=16; %number of grid points in xyz-directions
h=L/Nx;
wsearch=5; %Interpolation num, 5-grids for each direction.
width=wsearch*h; %Interpolation distance
width_window=width+0.5*h;
beta=25; %5 10 15 20
rx(1:Nx+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Nx+2*wsearch+2); %With virtual images
ry(1:Ny+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Ny+2*wsearch+2); %With virtual images
rz(1:Nz+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Nz+2*wsearch+2);

pad_ratio=1000; %Padding ratio w.r.t window width
pad_search=pad_ratio*wsearch;
pad_width=pad_search*h;
rz_padding(1:Nz+2*pad_search+2)=linspace(-L/2-h-pad_width,L/2+pad_width,Nz+2*pad_search+2);




%Precalculate scaling factor
TdK(1:Nx,1:Ny,1:2*pad_search+2)=0;
for i=1:Nx
    for j=1:Ny
        for k=1:Nz+2*pad_search+2
             fx=i-1; fy=j-1; fz=k-1;
            if (fx>ceil(Nx/2))
                fx=fx-Nx;
            end
            if (fy>ceil(Ny/2))
                fy=fy-Ny;
            end
            if (fz>(ceil(Nz/2)+pad_search+1))
                fz=fz-Nz-2*pad_search-2;
            end    
            k_x=fx*2*pi/L;
            k_y=fy*2*pi/L;
            k_z=fz*2*pi/(L+2*pad_width+2*h);
            TdK(i,j,k)=sqrt(pi)/4*wl*sl^3*exp(-1*sl^2*(k_x^2+k_y^2+k_z^2)/4);
        end
    end
end




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

HH=sum(sum(sum(H_temp)));


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


H(1:Nx,1:Ny,1:Nz+2*pad_search+2)=0;
H(:,:,(pad_search-wsearch)+1:Nz+(pad_search+wsearch)+2)=H_temp(wsearch+2:Nx+wsearch+1,wsearch+2:Ny+wsearch+1,1:Nz+2*wsearch+2);



%% Step2: HF_RtK. Use FFT to get H_hat
H_kspace(1:Nx,1:Ny,1:Nz+2*pad_search+2)=fftn(H);

%% Step3: ScalingH. Take the scaling of eack k-mode
%TdK(1:Nx+2*wsearch,1:Ny+2*wsearch,1:Nz+2*wsearch)=0;
H_tilde_kspace(1:Nx,1:Ny,1:Nz+2*pad_search+2)=0;
% Calculating scaling factor
for i=1:Nx
    for j=1:Ny
        for k=1:Nz+2*pad_search+2
             fx=i-1; fy=j-1; fz=k-1;
            if (fx>ceil(Nx/2))
                fx=fx-Nx;
            end
            if (fy>ceil(Ny/2))
                fy=fy-Ny;
            end
            if (fz>(ceil(Nz/2)+pad_search+1))
                fz=fz-Nz-2*pad_search-2;
            end    
            k_x=fx*2*pi/L;
            k_y=fy*2*pi/L;
            k_z=fz*2*pi/(L+2*pad_width+2*h);
            H_tilde_kspace(i,j,k)=TdK(i,j,k)*(FWkb(k_x,width_window,beta)*FWkb(k_y,width_window,beta)*FWkb(k_z,width_window,beta))^(-2)*H_kspace(i,j,k);
        end
    end
end


%% Step4:HF_KtR. Use IFFT to get H_tilde
H_tilde(1:Nx,1:Ny,1:Nz+2*pad_search+2)=ifftn(H_tilde_kspace);

%% Step5:Numerical Integral for potential
phi(1:N3)=0; %Potential
for j=1:N3 %Contribute the j-th charge
    xj=x(j,1); yj=x(j,2); zj=x(j,3);
    
    % Find the neareast grid num
    xnear=floor((xj+L/2+0.5*h)/h)+1;
    ynear=floor((yj+L/2+0.5*h)/h)+1;
    znear=floor((zj+L/2+h+pad_width+0.5*h)/h)+1;
    
     for x_ind=-1*wsearch:wsearch
        for y_ind=-1*wsearch:wsearch
            for z_ind=-1*wsearch:wsearch
                
                z_interp=z_ind+znear;
                z_dist=abs(rz_padding(z_interp)-zj);
                
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

t_mid_padding=toc;
%% Verification
%phi
tic
realphi(1:N3)=0;
cutoff=30;

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

t_Dmid=toc;

%realphi
error=abs(phi-realphi);
acc=max(error);
% end
% acc