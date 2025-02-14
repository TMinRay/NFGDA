function [groups, areas, centroids] = cluster_points(x, y, points, radius, show_fig)
    function [area, positions_sum] = dfs(r,c,color)
        % Process current node
        groups(r,c) = color;
        points(r,c) = false;
        area = 1;
        positions_sum = [0;0];

        % Find neighborhood
        curr_X = x(r,c);
        curr_Y = y(r,c);
        positions_sum = positions_sum + [curr_X; curr_Y];

        [~,c_low] = findclosest(x_coords, curr_X - radius/2);
        [~,c_high] = findclosest(x_coords, curr_X + radius/2);
        [~,r_low] = findclosest(y_coords, curr_Y - radius/2);
        [~,r_high] = findclosest(y_coords, curr_Y + radius/2);

        % Process children
        for nr = r_low: r_high
            for nc = c_low: c_high
                newX = x(nr, nc);
                newY = y(nr, nc);
                if (newX - curr_X)^2 + (newY - curr_Y)^2 <= radius^2 && points(nr,nc) && groups(nr, nc) == 0
                    [new_area, new_positions_sum] = dfs(nr, nc, color);
                    area = area + new_area;
                    positions_sum = positions_sum + new_positions_sum;
                end
            end
        end
    end

[NROWS, NCOLS] = size(x);
x_coords = x(1,:);
y_coords = y(:,1);

groups = zeros(NROWS, NCOLS);
n_groups = 0;
areas = [];
centroids = [];
for i = 1:NROWS
    for j = 1:NCOLS
        if points(i,j)
            n_groups = n_groups + 1;
            [area, positions_sum] = dfs(i,j, n_groups);
            areas = [areas, area];
            centroids = [centroids, positions_sum/area];
        end
    end
end


% Visualize clusters
if show_fig
    fig = figure;
    set(fig,'Position',[100 100 500 480]);

    markers = {'o','+','*', 'diamond', 'square', "^", "v"};
    for id = 1: length(areas)
        mask = groups == id;
        centroid = centroids(:,id);
        scatter(x(mask)', y(mask)', 'Marker',markers{1 + mod(id-1,7)}); hold on;
        plot(centroid(1), centroid(2),"r*"); hold on;
    end
    hold off
    xlim([-100,100]);
    ylim([-100,100]);
    title(sprintf("Clusters. Radius = %4.1f km", radius));
    grid on;
    set(gca,'TickDir','out','box','on','TickLength'  , [.01 .01], 'LineWidth' , 2);
    % print("clusters","-dpng");
end

end

function [val, idx] = findclosest(arr, target)
[val, idx] =  min(abs(arr - target));
end

