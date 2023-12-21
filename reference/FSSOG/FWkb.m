function wk=FWkb(k,width,beta)
wk=2*width*sinh(sqrt(beta^2-k^2*width^2))/(besseli(0,beta)*sqrt(beta^2-k^2*width^2));