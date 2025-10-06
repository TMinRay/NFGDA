clear all
close all

fuzzGST=readfis('./NF00ref_YHWANG.fis');
plotInputMFs(fuzzGST)

casename = 'KAMA20240227_21'
load(['../mat/train/train_set_',casename,'.mat']);
idx = randperm(size(trainNF, 1));

data_shuffled = trainNF(idx, :);
[Xn, ps] = mapminmax(data_shuffled', 0, 1);   % X is N-by-6, transpose to 6-by-N
data_shuffled = Xn';
radius = 0.5;
initFIS = genfis2(data_shuffled(:,1:6), data_shuffled(:,7), radius);
figure
plotInputMFs(initFIS)

opt = anfisOptions('InitialFIS',initFIS,'EpochNumber',250);
opt.ValidationData = data_shuffled(floor(size(data_shuffled,1)/10*4)+1:end,:);
[fis, trainError, stepSize, chkFIS, chkError] = anfis( ...
    data_shuffled(1:floor(size(data_shuffled,1)/10*4),:),opt);
figure
plotInputMFs(fis)



Xn = mapminmax('apply', trainNF', ps)';

figure
plot(evalfis(fis,Xn(:,1:end-1)),'LineWidth', 2)
hold on
plot(evalfis(initFIS,Xn(:,1:end-1)),'LineWidth', 2)
plot(Xn(:,end),'LineWidth', 2)

opt = anfisOptions('InitialFIS',fis,'EpochNumber',1000);
opt.ValidationData = data_shuffled(floor(size(data_shuffled,1)/10*4)+1:end,:);
[fis, trainError, stepSize, chkFIS, chkError] = anfis( ...
    data_shuffled(1:floor(size(data_shuffled,1)/10*4),:),opt);