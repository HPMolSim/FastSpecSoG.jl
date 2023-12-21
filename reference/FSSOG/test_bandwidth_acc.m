
L=100;
sl(1:10)=0;
b=1.62976708826776469; %Base 
sigma=3.633717409009413; %bandwidth
for mm=0:9
sl(mm+1)=sqrt(2)*b^mm*sigma;
end

ratio=1/5;
criterion=L*ratio;

acc_mid=load('acc_mid.mat','acc');
a_mid=struct2cell(acc_mid);
aa_mid=cell2mat(a_mid);
acc_wide=load('acc_wide.mat','acc_wide_store');
a_wide=struct2cell(acc_wide);
aa_wide=cell2mat(a_wide);

plot(sl,aa_mid,'r');
hold on;
scatter(sl,aa_mid,'r','filled');
hold on;
plot(sl,aa_wide,'b');
hold on;
scatter(sl,aa_wide,'b','filled');
hold on;
plot([criterion,criterion],[1e-16,0.2],'g--');
hold on;

legend("3DFFT","Data of 3DFFT","FFCT","Data of FFCT","Criteria: sl=L*1/5");

