Nx=[8,16,32,64,128];
acc=[5.42e-1,1.5e-1,3.95e-4,1.77e-8,1.92e-9];

plot(Nx,acc,'r');
hold on;
scatter(Nx,acc,'r','filled');
hold on;

legend("3DFFT","Data of 3DFFT");