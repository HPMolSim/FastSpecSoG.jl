% function  acc=FSSOG_mid_padding_opt(M,L,Nx,Ny,Nz,wsearch,pad_ratio,realphi)
%%  Initialization

% acc(1:10)=0;
% for mm=0:9
%u-series setting——Single Gaussian term
M=3; %Test Gaussian term
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
wl=(2*log(b))/sqrt(2*pi*sigma^2)*(1/b^M);
sl=sqrt(2)*b^M*sigma;



%Particle information(z-direction is free):(x,y,z,q)
L=100;
N=8;
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
Nx=16; Ny=Nx; Nz=Nx; %number of grid points in xyz-directions
h=L/Nx;
wsearch=Nx/4; % 窗函数宽度
width=wsearch*h; %Interpolation distance
width_window=width+0.5*h;

beta=2*2.5*wsearch; % kB function 参数
rx(1:Nx+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Nx+2*wsearch+2); %With virtual images
ry(1:Ny+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Ny+2*wsearch+2); %With virtual images
% rz(1:Nz+2*wsearch+2)=linspace(-L/2-h-width,L/2+width,Nz+2*wsearch+2);

pad_ratio=3; %Padding ratio w.r.t window width
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

kx(1:Nx)=linspace(0,2*pi*(Nx-1)/L,Nx);
ky(1:Ny)=linspace(0,2*pi*(Ny-1)/L,Ny);
kz(1:Nz+2*pad_search+2)=linspace(0,2*pi*(Nz+2*pad_search+2-1)/(L+2*pad_width+2*h),Nz+2*pad_search+2);
for i=1:Nx
    if (i-1>ceil(Nx/2))
        kx(i)=kx(i)-2*pi*Nx/L;
    end
end
for j=1:Ny
    if (j-1>ceil(Ny/2))
        ky(j)=ky(j)-2*pi*Ny/L;
    end
end
for k=1:Nz+2*pad_search+2
   if (k-1>(ceil(Nz/2)+pad_search+1))
       kz(k)=kz(k)-2*pi*(Nz+2*pad_search+2)/(L+2*pad_width+2*h);
   end  
end

FWKB_x(1:Nx)=0;
FWKB_y(1:Ny)=0;
FWKB_z(1:Nz+2*pad_search+2)=0;
for i=1:Nx
    FWKB_x(i)=FWkb(kx(i),width_window,beta);
end
for j=1:Ny
    FWKB_y(j)=FWkb(ky(j),width_window,beta);
end
for k=1:Nz+2*pad_search+2
    FWKB_z(k)=FWkb(kz(k),width_window,beta);
end


%Precompute Chebyshev polynomials
nu=10;
Coef(1:2*wsearch+2,1:nu)=PKB_coef(h,2*wsearch,width,beta,nu);


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

    %Store the interpolation index
    xx_interp(1:2*wsearch+1)=0;
    yy_interp(1:2*wsearch+1)=0;
    zz_interp(1:2*wsearch+1)=0;
    for ind=-1*wsearch:wsearch
        xx_interp(ind+wsearch+1)=ind+xnear;
        yy_interp(ind+wsearch+1)=ind+ynear;
        zz_interp(ind+wsearch+1)=ind+znear;
    end
    
    %Store the distance
    xx_dist(1:2*wsearch+1)=0;
    yy_dist(1:2*wsearch+1)=0;
    zz_dist(1:2*wsearch+1)=0;
    for ind=-1*wsearch:wsearch
        xx_dist(ind+wsearch+1)=abs(rx(xx_interp(ind+wsearch+1))-xj);
        yy_dist(ind+wsearch+1)=abs(ry(yy_interp(ind+wsearch+1))-yj);
        zz_dist(ind+wsearch+1)=abs(rz(zz_interp(ind+wsearch+1))-zj);
    end
    
    %Store the coeeficients
    Cx_int(1:2*wsearch+1,1:nu)=0;
    Cy_int(1:2*wsearch+1,1:nu)=0;
    Cz_int(1:2*wsearch+1,1:nu)=0;
    index_xx(1:2*wsearch+1)=0;
    index_yy(1:2*wsearch+1)=0;
    index_zz(1:2*wsearch+1)=0;
    for ind=-1*wsearch:wsearch
        index_xx(ind+wsearch+1)=floor((xx_dist(ind+wsearch+1)-(-width))/h)+2;
        index_yy(ind+wsearch+1)=floor((yy_dist(ind+wsearch+1)-(-width))/h)+2;
        index_zz(ind+wsearch+1)=floor((zz_dist(ind+wsearch+1)-(-width))/h)+2;
        if (abs(xx_dist(ind+wsearch+1))<=width_window) 
           Cx_int(ind+wsearch+1,:)=Coef(index_xx(ind+wsearch+1),:);
        end
        
        if (abs(yy_dist(ind+wsearch+1))<=width_window) 
           Cy_int(ind+wsearch+1,:)=Coef(index_yy(ind+wsearch+1),:);
        end
        
        if (abs(zz_dist(ind+wsearch+1))<=width_window) 
           Cz_int(ind+wsearch+1,:)=Coef(index_zz(ind+wsearch+1),:);
        end
    end
    
    %Calculate factor
    PKB_x(1:2*wsearch+1)=0;
    PKB_y(1:2*wsearch+1)=0;
    PKB_z(1:2*wsearch+1)=0;
    for ind=-1*wsearch:wsearch
        PKB_x(ind+wsearch+1)=PKB_calc(xx_dist(ind+wsearch+1),width,Cx_int(ind+wsearch+1,:),index_xx(ind+wsearch+1),2*wsearch,h,nu);
        PKB_y(ind+wsearch+1)=PKB_calc(yy_dist(ind+wsearch+1),width,Cy_int(ind+wsearch+1,:),index_yy(ind+wsearch+1),2*wsearch,h,nu);
        PKB_z(ind+wsearch+1)=PKB_calc(zz_dist(ind+wsearch+1),width,Cz_int(ind+wsearch+1,:),index_zz(ind+wsearch+1),2*wsearch,h,nu); 
    end   
    
    %Interpolate charge
    for x_ind=-1*wsearch:wsearch
        for y_ind=-1*wsearch:wsearch
            for z_ind=-1*wsearch:wsearch
                if (abs(xx_dist(x_ind+wsearch+1))>width_window) || (abs(yy_dist(y_ind+wsearch+1))>width_window) || (abs(zz_dist(z_ind+wsearch+1))>width_window)
                    continue;
                end
                H_temp(xx_interp(x_ind+wsearch+1),yy_interp(y_ind+wsearch+1),zz_interp(z_ind+wsearch+1))=H_temp(xx_interp(x_ind+wsearch+1),yy_interp(y_ind+wsearch+1),zz_interp(z_ind+wsearch+1))+x(j,4)*PKB_x(x_ind+wsearch+1)*PKB_y(y_ind+wsearch+1)*PKB_z(z_ind+wsearch+1);
            end
        end
    end
end

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
toc


%% Step2: HF_RtK. Use FFT to get H_hat
tic
H_kspace(1:Nx,1:Ny,1:Nz+2*pad_search+2)=fftn(H);
toc

%% Step3: ScalingH. Take the scaling of eack k-mode
tic
%TdK(1:Nx+2*wsearch,1:Ny+2*wsearch,1:Nz+2*wsearch)=0;
H_tilde_kspace(1:Nx,1:Ny,1:Nz+2*pad_search+2)=0;
% Calculating scaling factor

for i=1:Nx
    for j=1:Ny
        for k=1:Nz+2*pad_search+2
           H_tilde_kspace(i,j,k)=TdK(i,j,k)*(FWKB_x(i)*FWKB_y(j)*FWKB_z(k))^(-2)*H_kspace(i,j,k);
        end
    end
end
toc

%% Step4:HF_KtR. Use IFFT to get H_tilde
tic
H_tilde(1:Nx,1:Ny,1:Nz+2*pad_search+2)=ifftn(H_tilde_kspace);
toc

%% Step5:Numerical Integral for potential
tic
phi(1:N3)=0; %Potential
for j=1:N3 %Contribute the j-th charge
    xj=x(j,1); yj=x(j,2); zj=x(j,3);
    
    % Find the neareast grid num
    xnear=floor((xj+L/2+0.5*h)/h)+1;
    ynear=floor((yj+L/2+0.5*h)/h)+1;
    znear=floor((zj+L/2+h+pad_width+0.5*h)/h)+1;
    
    %Store the interpolation index
    xx_interp(1:2*wsearch+1)=0;
    yy_interp(1:2*wsearch+1)=0;
    zz_interp(1:2*wsearch+1)=0;
    xxx_temp(1:2*wsearch+1)=0;
    yyy_temp(1:2*wsearch+1)=0;
    for ind=-1*wsearch:wsearch
        xxx_temp(ind+wsearch+1)=ind+xnear;
        xx_interp(ind+wsearch+1)=mod(ind+xnear,Nx);
        if (xx_interp(ind+wsearch+1)==0)
             xx_interp(ind+wsearch+1)=Nx;
        end
        
        yyy_temp(ind+wsearch+1)=ind+ynear;
        yy_interp(ind+wsearch+1)=mod(ind+ynear,Ny);
        if (yy_interp(ind+wsearch+1)==0)
             yy_interp(ind+wsearch+1)=Ny;
        end

        zz_interp(ind+wsearch+1)=ind+znear;
    end
    
    %Store the distance
    xx_dist(1:2*wsearch+1)=0;
    yy_dist(1:2*wsearch+1)=0;
    zz_dist(1:2*wsearch+1)=0;
    for ind=-1*wsearch:wsearch
        xx_dist(ind+wsearch+1)=abs((-L/2+(xxx_temp(ind+wsearch+1)-1)*h)-xj);
        yy_dist(ind+wsearch+1)=abs((-L/2+(yyy_temp(ind+wsearch+1)-1)*h)-yj);
        zz_dist(ind+wsearch+1)=abs(rz_padding(zz_interp(ind+wsearch+1))-zj);
    end
    
    %Store the coeeficients
    Cx_int(1:2*wsearch+1,1:nu)=0;
    Cy_int(1:2*wsearch+1,1:nu)=0;
    Cz_int(1:2*wsearch+1,1:nu)=0;
    index_xx(1:2*wsearch+1)=0;
    index_yy(1:2*wsearch+1)=0;
    index_zz(1:2*wsearch+1)=0;
    for ind=-1*wsearch:wsearch
        index_xx(ind+wsearch+1)=floor((xx_dist(ind+wsearch+1)-(-width))/h)+2;
        index_yy(ind+wsearch+1)=floor((yy_dist(ind+wsearch+1)-(-width))/h)+2;
        index_zz(ind+wsearch+1)=floor((zz_dist(ind+wsearch+1)-(-width))/h)+2;
        if (abs(xx_dist(ind+wsearch+1))<=width_window) 
           Cx_int(ind+wsearch+1,:)=Coef(index_xx(ind+wsearch+1),:);
        end
        
        if (abs(yy_dist(ind+wsearch+1))<=width_window) 
           Cy_int(ind+wsearch+1,:)=Coef(index_yy(ind+wsearch+1),:);
        end
        
        if (abs(zz_dist(ind+wsearch+1))<=width_window) 
           Cz_int(ind+wsearch+1,:)=Coef(index_zz(ind+wsearch+1),:);
        end
    end

    %Calculate factor
    PKB_x(1:2*wsearch+1)=0;
    PKB_y(1:2*wsearch+1)=0;
    PKB_z(1:2*wsearch+1)=0;
    for ind=-1*wsearch:wsearch
        PKB_x(ind+wsearch+1)=PKB_calc(xx_dist(ind+wsearch+1),width,Cx_int(ind+wsearch+1,:),index_xx(ind+wsearch+1),2*wsearch,h,nu);
        PKB_y(ind+wsearch+1)=PKB_calc(yy_dist(ind+wsearch+1),width,Cy_int(ind+wsearch+1,:),index_yy(ind+wsearch+1),2*wsearch,h,nu);
        PKB_z(ind+wsearch+1)=PKB_calc(zz_dist(ind+wsearch+1),width,Cz_int(ind+wsearch+1,:),index_zz(ind+wsearch+1),2*wsearch,h,nu); 
    end   
    
    for x_ind=-wsearch:wsearch
        for y_ind=-wsearch:wsearch
            for z_ind=-wsearch:wsearch
                if (abs(xx_dist(x_ind+wsearch+1))>width_window) || (abs(yy_dist(y_ind+wsearch+1))>width_window) || (abs(zz_dist(z_ind+wsearch+1))>width_window)
                    continue;
                end
                phi(j)=phi(j)+4*pi*h^3*H_tilde(xx_interp(x_ind+wsearch+1),yy_interp(y_ind+wsearch+1),zz_interp(z_ind+wsearch+1))*PKB_x(x_ind+wsearch+1)*PKB_y(y_ind+wsearch+1)*PKB_z(z_ind+wsearch+1);              
            end
        end
    end
   
end
toc

tic
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