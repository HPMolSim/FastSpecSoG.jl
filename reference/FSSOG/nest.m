%Horner method to calculate polynomial
function y=nest(d,c,x)
y=c(d+1);
for i=d:-1:1
    y=y.*x+c(i);
end