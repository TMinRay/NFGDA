% TODO(pjatau) Add ground truth history
function [pos_history, ebox_history] = initialize_tracks(init_pos, init_vel, ts)
last_pos = init_pos;
vel_cart = [init_vel(1,:).*sin(init_vel(2,:)); init_vel(1,:).*cos(init_vel(2,:))];
last_t = 0;
pos_history{1} = init_pos;
ebox_history{1} = [];
for i = 2:length(ts)
    del_t = ts(i)- ts(i-1);
    del_pos = del_t*vel_cart;
    new_pos = last_pos + del_pos;

    %     noise = [rand(2,num_pts)*max(ts)*speed];
    noise = [];

    %     if rand(1) > 0
    %         clutter = [5+rand(1);5 + rand(1)];
    %     else
    %         clutter = [];
    %     end
    clutter = [];

    pos_history{i} = [new_pos noise clutter];
    ebox_history{i} = [];
    last_pos = new_pos;
end
end