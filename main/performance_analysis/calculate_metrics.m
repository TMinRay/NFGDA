
clear; clc; close all;

figure_folder = './figures';
line_out_folder = './gf_and_box_lines';
stats_folder = './stats';

if ~isfolder(stats_folder)
    mkdir(stats_folder);
end

color_order = get(gca,'colororder'); close;

case_names = {'KABX20210708_23'};
% case_names = {'KABX20210708_00_cured_gc','KABX20210708_23_cured_gc','KABX20210707_00_cured_gc'};

plds = []; precs = [];

for case_id = 1:length(case_names)

    case_name = case_names{case_id};

    if case_name == "."|| case_name == ".."
        continue
    end

    fig_path = fullfile(figure_folder, case_name);
    if ~isfolder(fig_path)
        mkdir(fig_path);
    end

    % ppi_folder = fullfile(line_out_folder, case_name);
    % idx_in_box = dir(ppi_folder);
    % ppi_names = {idx_in_box.name};
    v6m_path = fullfile('../../python/tracking_points/nf_preds',case_name);
    v6m_list = {dir(fullfile(v6m_path,'nf_pred*_V06.mat')).name};

    stats_path = fullfile(stats_folder, [case_name '.mat']);
   
    case_plds = [];
    frame_nums_plds = [];
    case_precisions = [];
    frame_nums_precs = [];
    for i = 1:length(v6m_list)

        % ppi_name_ext = ppi_names{i};
        % if ppi_name_ext == "." || ppi_name_ext == ".."
        %     continue
        % end

        % obj = strsplit(ppi_name_ext,".");
        % ppi_num = obj{1};
        obj = strsplit(v6m_list{i},"_");
        ppi_num = obj{3};

        % load(fullfile(line_out_folder, case_name,ppi_name_ext));
        line_out_path = fullfile(line_out_folder, case_name);
        load(fullfile(line_out_path,['gf_lines' v6m_list{i}(8:end)]));
        [num_boxes, num_areas] = size(gf_line_hit);

        % Calculate PLD
        gf_pld = cell([num_boxes 1]);
        avg_pld = 0;
        hit_l_ppi = 0;
        for i_box = 1:num_boxes
            % hit length
            hit_l_box = 0;
            for id = 1:num_areas
                hit = gf_line_hit{i_box, id};
                if hit
                    hit_l = arclength(hit(1,:), hit(2,:)); % arclength(x,y)
                    hit_l_box = hit_l_box + hit_l;
                end
            end

            % box length
            box_ma = box_major_axis{i_box};
            box_l = arclength(box_ma(1,:), box_ma(2,:));

            % pld
            gf_pld{i_box} = hit_l_box/box_l;
            avg_pld = avg_pld + gf_pld{i_box}/num_boxes;

            hit_l_ppi = hit_l_ppi + hit_l_box;
        end
        case_plds = [case_plds avg_pld];
        frame_nums_plds = [frame_nums_plds str2num(ppi_num)];

        % Precision
        total_pred_l = 0;
        for id = 1:num_areas
            pred_line = gf_lines{id};
            if isempty(pred_line)
                continue
            end

            pred_l = arclength(pred_line(1,:), pred_line(2,:));
            total_pred_l = total_pred_l + pred_l;
        end
        precision_ppi = hit_l_ppi/total_pred_l;
        case_precisions = [case_precisions, precision_ppi];
        frame_nums_precs = [frame_nums_precs str2num(ppi_num)];

        %% Figures
        fig = figure;
        set(fig,'Position',[100 100 500 480]);

        for id = 1:num_areas
            curr_points = gf_points{id};
            curr_line = gf_lines{id};

            if isempty(curr_points)
                continue
            end

            curr_color = color_order(1 + mod(id-1,length(color_order)),:);

            scatter(curr_points(1,:),curr_points(2,:),80,curr_color, "filled","AlphaData",1); hold on;
            alpha(0.3); hold on;
            plot(curr_line(1,:), curr_line(2,:), "LineWidth",3,"Color","r"); hold on;

            % hits
            for i_box = 1:num_boxes
                curr_line_hit = gf_line_hit{i_box, id};
                if curr_line_hit
                    plot(curr_line_hit(1,:), curr_line_hit(2,:),"LineWidth",3,"Color","k"); hold on;
                end
            end
        end

        % plot major axes of evaluation boxes
        contour(xi2,yi2,evalbox,'Color',"#EDB120",'linewidth',1); hold on;

        for i_box = 1:num_boxes
            box_ma = box_major_axis{i_box};

            %             if box_ma
            %                 plot(box_ma(1,:) , box_ma(2,:) ,"LineWidth",3,"Color","m"); hold on;
            %             end
            text(box_ma(1,1), box_ma(2,1),sprintf("%4.2f", gf_pld{i_box}),FontSize=14);
        end

        xlim([-90,90]);
        ylim([-90,90]);
        grid on;
        title(sprintf("Average PLD: %4.2f. Precision: %4.2f.",avg_pld, precision_ppi),'Fontsize',14);
        set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth' , 2);
        % print(fullfile(fig_path,[ppi_num '.png']),'-dpng')
        print(fullfile(fig_path,['gf_lines' v6m_list{i}(8:end-4) '.png']),'-dpng')

    end
    fprintf("Averages for %s. PLD: %4.2f. Precision: %4.2f\n", case_name, nanmean(case_plds), nanmean(case_precisions));

    plds = [plds case_plds];
    precs = [precs case_precisions];
    save(stats_path,'case_plds','case_precisions','frame_nums_plds','frame_nums_precs');
end

nanmean(plds)
nanmean(precs)

