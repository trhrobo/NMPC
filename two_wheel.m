goal_error = 1.0e-2;
%[x y theta]
now_pos = [0 0 0];
goal_pos = [100 100 0];

dt=0.01;
iteration_time=30;
iteration_num=iteration_time/dt;
%現在の状態
%x={x1, x2}
init_X=[-4.5; 1.5; 0.25];
disp("pp")
nmpc = NMPC(init_X);
%X(0)を測定する(初期値を代入する)
for i = 1:iteration_num-1
    time=i*dt;
    u=nmpc.CGMRES(time);
    nmpc.updateState(u, dt);
end
plot(nmpc.save_x(:,1), nmpc.save_x(:,2))

