function val=SOG(Bcut,w0,base,sigma,r)

% w0=sqrt(2*pi)*sigma/(2*log(base)*r0*exp(-r0^2/(2*sigma^2)));
% for i=1:Bcut
%     w0=w0-(1/exp(-r0^2/(2*sigma^2)))*(1/base^i)*(exp(-(r0/(base^i*sigma))^2/2));
% end

temp=2*log(base)*w0/(sqrt(2*pi)*sigma)*exp(-r.^2/(2*sigma^2));
for i=1:Bcut
    temp=temp+2*log(base)/(base^i*sqrt(2*pi*sigma^2))*exp(-r.^2/(2*base^(2*i)*sigma^2));
end

val=temp;