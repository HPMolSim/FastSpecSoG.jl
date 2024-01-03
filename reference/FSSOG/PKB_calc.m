function wr=PKB_calc(r,width,C,index,P,h,nu)

if (index>=2) && (index<=P+1)
    center=-width+(index-3/2)*h;
    scale=1/2*h;
    r_calc=(r-center)/scale;
end

if (index==1)
    center=-width-1/4*h;
    scale=1/4*h;
    r_calc=(r-center)/scale; 
end

if (index==P+2)
    center=width+1/4*h;
    scale=1/4*h;
    r_calc=(r-center)/scale; 
end

wr=nest(nu-1,C,r_calc);
