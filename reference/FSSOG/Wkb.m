function wr=Wkb(r,width,beta)
% if (abs(r)<=width)
%     wr=besseli(0,beta*sqrt(1-(r/width)^2))/besseli(0,beta);
% else
%     wr=0;
% end
wr=besseli(0,beta.*sqrt(1-(r./width).^2))./besseli(0,beta).*(abs(r)<=width);


