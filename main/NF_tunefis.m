clear all
close all

fuzzGST=readfis('./NF00ref_YHWANG.fis');
plotInputMFs(fuzzGST)

casename = 'KABX20230720_23'
load(['../mat/train/train_set_',casename,'.mat']);
idx = randperm(size(trainNF, 1));

data_shuffled = trainNF(idx, :);
% opt = anfisOptions('InitialFIS',fuzzGST,'EpochNumber',20);
% fis = anfis(data_shuffled,opt);
% [fis,trainFISError] = anfis(data_shuffled,opt);
fises={};
errors={};
fis=fuzzGST;
% for i=1:5
%     opt = anfisOptions('InitialFIS',fis,'EpochNumber',2000);
%     [fis,trainFISError] = anfis(data_shuffled,opt);
%     fises{i}=fis;
%     errors{i}=trainFISError;
% end
opt = anfisOptions('InitialFIS',fis,'EpochNumber',200);
opt.ValidationData = data_shuffled(floor(size(data_shuffled,1)/10*4)+1:end,:);
[fis, trainError, stepSize, chkFIS, chkError] = anfis( ...
    data_shuffled(1:floor(size(data_shuffled,1)/10*4),:),opt);

figure
plotInputMFs(fis)

figure
x = [1:size(trainError,1)];
plot(x, trainError,'.b',x, chkError,'*r');
xlabel('Epochs (no units)')
ylabel('Loss (unknown)')
legend train val
grid on
% fis.Name = ['tune ' casename];
% writeFIS(fis, './tuned.fis');
% 
% set(groot, 'DefaultLineLineWidth', 2);        % Thicker lines
% set(groot, 'DefaultAxesFontSize', 14);        % Larger axis labels
% set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.2);  % Scale label font
% set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.3);  % Scale title font
% set(groot, 'DefaultLineMarkerSize', 8);  
% 
% hold on
% for i=1:5
% plot(2000*(i-1)+1:2000*i,errors{i})
% end
% ylabel('Error', 'FontSize', 14)
% xlabel('epoch', 'FontSize', 14)
% 
% figure
% hold on
% for i=1:5
%     plot(evalfis(fises{i},trainNF(:,1:end-1)),'LineWidth', 2)
% end
% warning('off','fuzzy:general:diagEvalfis_OutOfRangeInput')  % turn off all warnings (not recommended generally)
% plot(evalfis(fuzzGST,trainNF(:,1:end-1)),'LineWidth', 2)
% warning('on','fuzzy:general:diagEvalfis_OutOfRangeInput')
% plot(trainNF(:,end),'LineWidth', 2)
% legend('2000','4000','6000','8000','10000','no tune','truth')
% ylabel('predict', 'FontSize', 14)
% xlabel('sample', 'FontSize', 14)

figure
labels = {'beta', 'Z', 'rhoHV', ...
          'ZDR', 'SDv', 'SDphi'};
for iv=1:6
    subplot(2,3,iv); % 2 rows, 3 columns, index from 1 to 6
    histogram(trainNF(trainNF(:,7)==0,iv), 30, 'FaceColor','blue', 'FaceAlpha',0.5);
    hold on
    histogram(trainNF(trainNF(:,7)==1,iv), 30, 'FaceColor','red',  'FaceAlpha',0.5);
    title(labels{iv});
end
