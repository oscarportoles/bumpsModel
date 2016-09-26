
function [lkhst2] = generateDataDiscoverDiffStates(SNR,normedscore10,prefix10)
%% INPUT and OUTPUT
% normedscore10 - EEG data on each of the 10 PCAs
% prefix10 - index for subject(1),condition(3)
% lkhst2 - likelihood under different number of states
%% main
num = 8;
lkhs = zeros(num);
mags = cell(num);
params = cell(num);
signal = cell(num);
prefix = cell(num);
x2=find(prefix10(:,3)==1);

lkhst1 = zeros(num,8,20);
lkhst2 = zeros(num,8,20);
lkhst3 = zeros(num,8,20);
magst = cell(num,8,20);
paramst = cell(num,8,20);
eventprobst = cell(num,8,20);
for i = 1:num % number of actual states

    [lkhs(i) mags{i} params{i} eventprobs{i}]=hsmmEEG(normedscore10,zeros(10,i),repmat([1 50],i+1,1),1,prefix10,2,0);

    %[lkhN magsN paramsN signal prefix] = runnoisenew(size(x2,1),mags{i},params{i},SNR);
    [signal prefix pos pure LP] = genNoiseRegular(size(x2,1),mags{i},params{i}, SNR);
    x=find(prefix(:,3)==1);
    N = size(x,1); 
    y = x-1;
    y(1) = [];
    y = [y;size(prefix,1)];
    lens = y-x+1;
    chunk = floor(N/20);

    for h = 1:8 % number of states to look at
        for k = 1:1 % 1 fold is representive of 20 folds in LOOCV for synthetic data (sysmetric)
             
             tmp1 = h
             tmp2 = k
             a = x(chunk*(k-1)+1);
             b = y(chunk*k);
             if(k==20) 
                 b = y(end); 
             end
             
             data_train = signal;
             data_train(a:b,:) = [];
             data_test = signal(a:b,:);
             
             ind_train = prefix;
             ind_train(a:b,:) = [];
             ind_test = prefix(a:b,:);

             [lkhst1(i,h,k) magst{i,h,k} paramst{i,h,k} eventprobst{i,h,k}]=hsmmEEG(data_train,zeros(10,h),repmat([1 50],h+1,1),1,ind_train,2,0); 
             [lkhst2(i,h,k) magst{i,h,k} paramst{i,h,k} eventprobst{i,h,k}]=hsmmEEG(data_test,magst{i,h,k},paramst{i,h,k},0,ind_test,2,0);
             [lkhst3(i,h,k) magst{i,h,k} paramst{i,h,k} eventprobst{i,h,k}]=hsmmEEG(data_test,zeros(10,h),paramst{i,h,k},0,ind_test,2,0);
        end
           

   end 
end


figure,

for i = 1:8
   subplot(2,4,i)
   plot(squeeze(lkhst2(i,:,1)))
   title([num2str(i),' peak'])
    
end


    
end

