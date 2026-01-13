% labels = {'beta', 'Z', 'rhoHV', ...
%           'ZDR', 'SDv', 'SDphi'};
% 
% for iv=1:6
%     subplot(2,3,iv); % 2 rows, 3 columns, index from 1 to 6
%     x = linspace(fuzzGST.Inputs(iv).range(1),fuzzGST.Inputs(iv).range(2), 100);
%     hold on
%     for im=1:2
%         c = fuzzGST.Inputs(iv).MembershipFunctions(im).Parameters(2);
%         sig = fuzzGST.Inputs(iv).MembershipFunctions(im).Parameters(1);
%         mu = exp(-0.5 * ((x - c) / sig).^2);
%         plot(x, mu, 'LineWidth', 2);
%     end
%     grid on;
%     % title(sprintf('var %d', iv));
%     title(labels{iv});
% end
% 
% 
% 
% % Gaussian membership function formula
% mu = exp(-0.5 * ((x - c) / sigma).^2);
% 
% plot(x, mu, 'LineWidth', 2);
% 
% figure
% plotInputMFs(fuzzGST)

function plotInputMFs(fis)
% plotInputMFs  Plot membership functions for each input of a FIS
%
%   plotInputMFs(fis) plots the first two membership functions
%   for each input variable of the fuzzy system fis in a 2x3 subplot.

    labels = {'beta', 'Z', 'rhoHV', ...
              'ZDR', 'SDv', 'SDphi'};

    nInputs = numel(fis.Inputs);

    for iv = 1:nInputs
        subplot(2, 3, iv); % 2 rows, 3 columns
        x = linspace(fis.Inputs(iv).Range(1), fis.Inputs(iv).Range(2), 100);
        hold on
        for im = 1:min(2, numel(fis.Inputs(iv).MembershipFunctions))
            c   = fis.Inputs(iv).MembershipFunctions(im).Parameters(2);
            sig = fis.Inputs(iv).MembershipFunctions(im).Parameters(1);
            mu  = exp(-0.5 * ((x - c) / sig).^2);
            plot(x, mu, 'LineWidth', 2);
        end
        grid on
        title(labels{iv}, 'FontWeight', 'bold')
    end
end