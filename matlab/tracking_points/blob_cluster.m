function [groups, areas, centroids] = cluster_points(x, y, points, radius, show_fig)
    function [] = dfs(p,color)
        % Process current node
        in_this = (((xp - xp(p)).^2 + (yp - yp(p)).^2) < radius^2) & gp==0;
        gp(in_this) = color;
        for bp = find(in_this).'
            dfs(bp,color)
        end
    end
npoints = sum(points,'all');
xp = x(points);
yp = y(points);
gp = zeros(npoints,1);
n_groups = 0;
for p = 1:npoints
    if gp(p)==0
        n_groups = n_groups + 1;
        dfs(p,n_groups)
    end
end

groups = zeros(size(x));
groups(points) = gp;

areas=zeros(1,n_groups);
centroids=zeros(2,n_groups);
mask = groups~=0;
xp = x(mask);
yp = y(mask);
gp = groups(mask);
for p = 1:length(gp)
    areas(gp(p)) = areas(gp(p)) + 1;
    centroids(:,gp(p)) = centroids(:,gp(p)) + [xp(p); yp(p)];
end
centroids = centroids./areas;

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