goal_error = 1.0e-2;
%[x y theta]
dt=0.01;
iteration_time=10;
iteration_num=iteration_time/dt;
%現在の状態
goal_pos=[40; 40; 0];
%x={x1, x2}
init_X=[0; 0; 0.25];
%init_X=[1.5; 1.5; 0.25];
disp("pp")

field = zeros(300, 300);
for x = 25:30
    for y = 25:30
        field(x, y) = 10;
    end
end
obs_pos = [10 10];
nmpc = NMPC_two_wheel_obs(init_X, goal_pos, field, obs_pos);
%X(0)を測定する(初期値を代入する)
for i = 1:iteration_num-1
    time=i*dt;
    disp(time)
    u=nmpc.CGMRES(time, goal_pos);
    nmpc.updateState(u, dt);
    nmpc.figGraph();
end
[~, curvature_nmpc, ~] = curvature(nmpc.save_x(:,1:2));
curvature_nmpc=1/curvature_nmpc;

tiledlayout(2, 1)

% Tile 1
nexttile
plot(nmpc.save_x(:,1), nmpc.save_x(:,2))
title("pos")
