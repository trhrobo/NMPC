goal_error = 1.0e-2;
%[x y theta]
dt=0.01;
iteration_time=30;
iteration_num=iteration_time/dt;
%現在の状態
goal_pos=[40; 40; 0];
%x={x1, x2}
init_X=[0; 0; 0.25];
%init_X=[1.5; 1.5; 0.25];
disp("pp")
nmpc = NMPC_two_wheel(init_X, goal_pos);
%X(0)を測定する(初期値を代入する)
for i = 1:iteration_num-1
    time=i*dt;
    disp(time)
    u=nmpc.CGMRES(time, goal_pos);
    nmpc.updateState(u, dt);
end
[~, curvature_nmpc, ~] = curvature(nmpc.save_x(:,1:2));
curvature_nmpc=1/curvature_nmpc;

tiledlayout(2, 1)

% Tile 1
nexttile
plot(nmpc.save_x(:,1), nmpc.save_x(:,2))
title("pos")

% Tile 2
nexttile
plot(curvature_nmpc)
title("curvature")