%field
x = [0;  0;  5; 10; 15; 17; 15; 20; 20; 20; 25; 30; 35; 40];
y = [0; 40;  5; 25; 10;  5; 0; 10; 25; 15;  5; 20;  8; 40];
z = [5; 10;  5; 10; 30; 30; 30;  5;  5;  5; 30;  5; 20; 10];
%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( x, y, z );

% Set up fittype and options.
ft = 'thinplateinterp';

% Fit model to data.
%[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

load fieldData.mat

%[x y theta]
dt=0.01;
iteration_time=10;
iteration_num=iteration_time/dt;
%現在の状態[x y theta]
goal_pos=[30; 20; 0];
%x={x1, x2}
init_X=[0; 0; 0.25];
%init_X=[1.5; 1.5; 0.25];
init_z = 5;
nmpc = NMPC3dem(init_X, goal_pos, init_z, resultMat);
%X(0)を測定する(初期値を代入する)
for i = 1:iteration_num-1
    time=i*dt;
    disp(time)
    u=nmpc.CGMRES(time, goal_pos, fitresult);
    nmpc.updateState(u, dt);
end
[~, curvature_nmpc, ~] = curvature(nmpc.save_x(:,1:2));
curvature_nmpc=1/curvature_nmpc;

tiledlayout(3, 1)


%Tile 1
nexttile
plot(fitresult)
title("field")

%Tile 2
nexttile
plot(fitresult)
hold on
plot3(nmpc.save_x(:,1), nmpc.save_x(:,2), nmpc.save_x(:,3), 'Color', 'k', 'LineWidth', 5)
hold off
title("pos")

%Tile 3
nexttile
plot(curvature_nmpc)
title("curvature")