L=100;
sl(1:10)=0;
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
for mm=1:10
sl(mm)=sqrt(2)*b^(mm+2)*sigma;
end

t_FFT(1:10)=[1.01 1.12 1.14 1.14 1.16 1.18 1.21 1.31 1.60 1.89];
t_FFCT(1:10)=1e-2.*[9.81 8.96 4.90 4.57 3.37 2.48 2.55 2.64 3.19 2.11];

plot(sl,t_FFT,'r');
hold on;
scatter(sl,t_FFT,'r','filled');
hold on;
plot(sl,t_FFCT,'b');
hold on;
scatter(sl,t_FFCT,'b','filled');
hold on;

legend("3DFFT","Data of 3DFFT","FFCT","Data of FFCT");
