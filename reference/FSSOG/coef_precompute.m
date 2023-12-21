function Coef=coef_precompute(delta,Lz,RP)

gamma(1:RP)=2;
gamma(1)=1;

Cnodes(1:RP)=0;
for i=1:RP
    Cnodes(i)=-Lz/2*cos((2*i-1)/(2*RP)*pi);
end

En(1:RP)=0;
for i=1:RP
    n=i-1;
    if (mod(n,2)==0)
      En(i)=gamma(i)*(1i)^(n/2)*exp(-1/(2*delta))*besselj(n/2,1i/(2*delta));
    end
end

a(1:RP,1:RP,1:RP)=0;
for k=1:RP
    for l=1:RP+1-k
        for n=k+l-1:RP
            for i=1:RP
                for j=1:RP
                    a(k,l,n)=a(k,l,n)+gamma(k)*gamma(l)/RP^2*ChebPoly(n-1,(Cnodes(i)-Cnodes(j))/2,Lz/2)*ChebPoly(k-1,Cnodes(i),Lz/2)*ChebPoly(l-1,Cnodes(j),Lz/2);
                end
            end
        end
    end
end

E(1:RP,1:RP)=0;
for k=1:RP
    for l=1:RP+1-k
        for n=k+l-1:RP
            E(k,l)=E(k,l)+a(k,l,n)*En(n);
        end
    end
end

Coef=E;