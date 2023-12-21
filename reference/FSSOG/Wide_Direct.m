function [real_phi,t_Dlong]=Wide_Direct(M_min,M_max,L,x,N)
M=M_min:1:M_max;
l_range=length(M);
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
wl=(2*log(b))./sqrt(2*pi*sigma^2).*(1./b.^M);
sl=sqrt(2).*b.^M*sigma;  
N3=N^3;

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

real_phi=realphi;
t_Dlong=toc;