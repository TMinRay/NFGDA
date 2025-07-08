clear; clc; close all;

data_path = './clusters';
figure_folder = './figures';
line_out_folder = './gf_and_box_lines';

color_order = get(gca,'colororder');
close;

case_names = {'KABX20210708_00_cured_gc'};
% case_names = {'KABX20210708_00_cured_gc','KABX20210708_23_cured_gc','KABX20210707_00_cured_gc'};

for case_id = 1:length(case_names)

    case_name = case_names{case_id};

    if case_name == "."|| case_name == ".."
        continue
    end

%     fig_path = fullfile(figure_folder, case_name);
%     if ~isfolder(fig_path)
%         mkdir(fig_path);
%     end

    line_out_path = fullfile(line_out_folder, case_name);
    if ~isfolder(line_out_path)
        mkdir(line_out_path);
    end

    ppi_folder = fullfile(data_path, case_name);
    idx_in_box = dir(ppi_folder);
    ppi_names = {idx_in_box.name};

    for i = 1:length(ppi_names)

        ppi_name_ext = ppi_names{i};
        if ppi_name_ext == "." || ppi_name_ext == ".."
            continue
        end

        obj = strsplit(ppi_name_ext,".");
        ppi_num = obj{1};

        ppi_path = fullfile(data_path, case_name, [ppi_num '.mat']);
        load(ppi_path);

        % Estimate GF lines
        gf_points = cell(size(areas));
        gf_lines = cell(size(areas));
        for id = 1:length(areas)
            mask = groups == id;
            if sum(sum(mask)) == 0
                continue
            end

            x = xi2(mask)';
            y = yi2(mask)';
            gf_points{id} = [x; y];

            [p, angle_rot, x_rot] = polyfit_rotation(x,y,2);
            curr_line = polyval_rotation(p, x_rot,angle_rot);
            gf_lines{id} = curr_line;
        end

        % Estimate major axes of evaluation box and intersection with GF line
        CC = bwconncomp(evalbox);
        num_boxes = CC.NumObjects;
        gf_line_hit = cell(num_boxes, length(areas));
        box_major_axis = cell([num_boxes 1]);

        for i_box = 1:num_boxes

            idx_box = CC.PixelIdxList{i_box};
            gfs_in_box = unique(groups(idx_box));

            x_box = xi2(idx_box)';
            y_box = yi2(idx_box)';

            % Obtain line accross eval box
            [p_box, angle_rot, x_rot] = polyfit_rotation(x_box,y_box,2);
            box_ma = polyval_rotation(p_box, x_rot,angle_rot);
            box_major_axis{i_box} = box_ma;

            for id = gfs_in_box'
                if id == 0
                    continue
                end

                curr_line = gf_lines{id};
                curr_line_quantized = round(curr_line./0.5)*0.5; % match to common grid

                % Find GF line segment that overlaps with eval box
                [idx_in_box]= inpolygon(curr_line_quantized(1,:),curr_line_quantized(2,:),x_box,y_box);
                line_in_box = curr_line(:,idx_in_box);
                gf_line_hit{i_box, id} = line_in_box;

            end
        end

        save(fullfile(line_out_path,ppi_name_ext),'xi2','yi2','evalbox','gf_points','gf_lines','box_major_axis','gf_line_hit');


        % Figures
        fig = figure;
        set(fig,'Position',[100 100 500 480]);

        for id = 1:length(areas)
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
        for i_box = 1:num_boxes
            box_ma = box_major_axis{i_box};
            if box_ma
                plot(box_ma(1,:), box_ma(2,:),"LineWidth",3,"Color","m"); hold on;
            end
        end

        contour(xi2,yi2,evalbox,'y-','linewidth',1); hold on;
        xlim([-80,80]);
        ylim([-80,80]);
        grid on;
        set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth' , 2);
%         print(fullfile(fig_path,[ppi_num '.png']),'-dpng')
    end
end


