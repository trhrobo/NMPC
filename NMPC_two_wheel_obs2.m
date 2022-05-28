classdef NMPC_two_wheel_obs2 < handle
    properties 
        %予測ステップ
        N_step=10;
        alpha=0.5;
        tf=1.0;
        zeta=100.0;
        ht=0.01;
        %x={x, y, theta}    
        x_size=3; 
        %各時刻における制御入力
        %・入力が2つ
        %u={u1, u2, dummy}
%         u_size=3;
        u_size=2;
        %U={u1, u2, dummy, mu}
%         U_size=4;
        U_size=2;
        mu;
        %f={rHru(x_size+dummy+mu)}
%         f_size=4;
        f_size=2;
        max_iteration;
        U;
        dt=0;
        save_x;
        save_u;
        goal_pos;
        X;
        X_cal;
        param_pos=1.0;
        param_pos_theta=1.0
        param_log=30.0;
        param_u_v=1;
        param_u_omega=0;
        param_dummy=0;
    end
    methods
        function obj = NMPC_two_wheel_obs2(X_, goal_pos_)
            obj.max_iteration = obj.U_size*obj.N_step;
            obj.U=zeros(obj.max_iteration, 1);
            obj.X=X_;
            obj.goal_pos=goal_pos_;
            for i = 0:obj.N_step-1
                obj.U((obj.U_size*i)+1, 1)=1.0;
                obj.U((obj.U_size*i)+2, 1)=0.5;
%                 obj.U((obj.U_size*i)+3, 1)=100-3^2;
%                 obj.U((obj.U_size*i)+4, 1)=0.011;
            end
        end
        function figGraph(obj)
            r=4;
            xc=0;
            yc=0;

            theta=linspace(0,2*pi);
            x=r*cos(theta)+xc;
            y=r*sin(theta)+yc;

            plot(x,y)

            X_=(reshape(obj.X_cal, obj.x_size, []))';
            %obj.figgraph = animatedline(obj.save_x(:,1), obj.save_x(:,2));
            figgraph1 = animatedline(X_(:,1), X_(:, 2));
            axis([-30, 50, -30, 50]);
            figgraph1.Color='green';
            figgraph2 = animatedline(obj.save_x(:, 1), obj.save_x(:, 2));
            axis([-30, 50, -30, 50]);
            figgraph2.Color='black';
            drawnow
        end
        function updateState(obj, U_, dt_)
            %状態Xを更新する
            %disp("dt_")
            %disp(dt_);
            %disp("dt_end")
            obj.X=obj.X+obj.calModel(obj.X, U_(1:obj.U_size,1))*dt_;
            %disp("U_")
            %disp(U_)
            %disp("U_end")
            disp([obj.X(1,1), obj.X(2,1)])
            obj.save_x=[obj.save_x;(obj.X)'];
            obj.save_u=[obj.save_u;(U_(1:obj.U_size,1))'];
        end
        function model_F = calModel(obj, X_, U_)
            %U={u, v, rho}
            model_F = [cos(X_(3, 1))*U_(1, 1);
                       sin(X_(3, 1))*U_(1, 1);
                       U_(2, 1)];
            %disp("model_F")
            %disp(model_F)
            %disp("model_F_end")
        end
        function rphirx = rphirx(obj, X_)
            x_=X_(1,1);
            y_=X_(2,1);
            theta_=X_(3,1);
            g_x_=obj.goal_pos(1,1);
            g_y_=obj.goal_pos(2,1);
            g_theta_=obj.goal_pos(3,1);
            %rphirx=[obj.param_pos*(x_-g_x_)-obj.param_log*2*x_/(x_*x_+y_*y_-16), obj.param_pos*(y_-g_y_)-obj.param_log*2*y_/(x_*x_+y_*y_-16), 0];
            rphirx=[obj.param_pos*(x_-g_x_)-obj.param_log*2*x_/(x_*x_+y_*y_-36), obj.param_pos*(y_-g_y_)-obj.param_log*2*y_/(x_*x_+y_*y_-36), obj.param_pos*(theta_-g_theta_)];
        end
        function rhrx = rHrx(obj, X_, u_, lamda_)
            x_=X_(1, 1);
            y_=X_(2, 1);
            theta_=X_(3, 1);
            u_v_=u_(1, 1);
            lamda1_=lamda_(1, 1);
            lamda2_=lamda_(2, 1);
            g_x_=obj.goal_pos(1,1);
            g_y_=obj.goal_pos(2,1);
            g_theta_=obj.goal_pos(3,1);
            %rhrx=[obj.param_pos*(x_-g_x_)-obj.param_log*2*x_/(x_*x_+y_*y_-16), obj.param_pos*(y_-g_y_)-obj.param_log*2*y_/(x_*x_+y_*y_-16),-obj.param_u_v*lamda1_*u_v_*sin(theta_)+obj.param_u_v*lamda2_*u_v_*cos(theta_)];
            rhrx=[obj.param_pos*(x_-g_x_)-obj.param_log*2*x_/(x_*x_+y_*y_-16), obj.param_pos*(y_-g_y_)-obj.param_log*2*y_/(x_*x_+y_*y_-16),obj.param_pos*(theta_-g_theta_)-obj.param_u_v*lamda1_*u_v_*sin(theta_)+obj.param_u_v*lamda2_*u_v_*cos(theta_)];
        end
        function cgmres = CGMRES(obj, time_, goal_pos_)
            obj.dt=obj.tf*(1-exp(-obj.alpha*time_))/obj.N_step;
            obj.goal_pos=goal_pos_;
            %gmres法を用いてdUを求める
            dU=zeros(obj.max_iteration, 1);
            gmres_Xm=zeros(obj.max_iteration, 1);
            gmres_X0=zeros(obj.max_iteration, 1);
            gmres_Ym=zeros(obj.max_iteration, 1);
            gmres_V=zeros(obj.max_iteration, obj.max_iteration+1);
            gmres_Vm=zeros(obj.max_iteration, obj.max_iteration);
            gmres_R0=zeros(obj.max_iteration, 1);
            %gmres_X0を設定する
            gmres_X0=obj.U;
            %disp("gmres_R0")
            %disp(gmres_R0)
            %disp("gmres_R0_end")
            %初期残差gmres_R0を求める
            gmres_R0=obj.calR0();
            gmres_V(:,1)=normc(gmres_R0);
            g=zeros(obj.max_iteration+1,1);
            g(1, 1)=norm(gmres_R0);
            h=zeros(obj.max_iteration+1, obj.max_iteration);
            Rm=zeros(obj.max_iteration, obj.max_iteration);
            c=zeros(obj.max_iteration, 1);
            s=zeros(obj.max_iteration, 1);
            temp_sigma=zeros(obj.U_size*obj.N_step, 1);
            tempAv=zeros(obj.U_size*obj.N_step, obj.max_iteration+1);
            tempcalAv=zeros(obj.U_size*obj.N_step, 1);
            for i = 1:obj.max_iteration
                temp_sigma=zeros(obj.U_size*obj.N_step, 1);
                tempcalAv=obj.calAv(gmres_V(:, i));
                for k = 1:i+1
                    tempAv(:, k)=tempcalAv;
                    h(k, i)=dot(tempAv(:, k), gmres_V(:, k));
                    temp_sigma=temp_sigma+h(k, i)*gmres_V(:, k);
                end
                temp_V=tempAv(:, i)-temp_sigma;
                h(i+1, i)=norm(temp_V);
                gmres_V(:, i+1)=temp_V/h(i+1, i);
                Rm(1, i)=h(1, i);
                for k = 1:i-1
                    temp1=c(k, 1)*Rm(k, i)+s(k, 1)*h(k+1, i);
                    temp2=-1*s(k, 1)*Rm(k, i)+c(k, 1)*h(k+1, i);
                    Rm(k, i)=temp1;
                    Rm(k+1, i)=temp2;
                end
                c(i, 1)=Rm(i, i)/sqrt(Rm(i, i)*Rm(i, i)+h(i+1, i)*h(i+1, i));
                s(i, 1)=h(i+1, i)/sqrt(Rm(i, i)*Rm(i, i)+h(i+1, i)*h(i+1, i));
                g(i+1, 1)=-1*s(i, 1)*g(i, 1);
                g(i, 1)=c(i, 1)*g(i, 1);
                Rm(i, i)=c(i, 1)*Rm(i, i)+s(i, 1)*h(i+1, i);
            end
            Gm=g(1:obj.max_iteration, 1);
            gmres_Vm=gmres_V(1:obj.max_iteration, 1:obj.max_iteration);
            %後退代入によってRm*Ym=Gmを解く
            for i = obj.max_iteration:-1:1
                temp_sigma_back=0;
                for k = i+1:obj.max_iteration
                    temp_sigma_back=temp_sigma_back+Rm(i, k)*gmres_Ym(k, 1);
                end
                gmres_Ym(i, 1)=(Gm(i, 1)-temp_sigma_back)/Rm(i, i);
            end
            gmres_Xm=gmres_X0+gmres_Vm*gmres_Ym;
            dU=gmres_Xm;
            obj.U=obj.U+dU*obj.ht;
            cgmres = obj.U(1:obj.U_size, 1);
        end      
        function F = calF(obj, U_, x_)
            F=zeros(obj.f_size*obj.N_step, 1);
            %制約なし
            %0~Nまでx1, x2
            X_=zeros(obj.x_size*(obj.N_step+1), 1);
            %1~Nまでlamda1, lamda2
            Lamda_=zeros(obj.x_size*obj.N_step, 1);
            %x_(x*)を求める
            %x0*(t)=x(t)を計算する
            X_(1:obj.x_size, 1)=x_;
            %xi+1*=xi*+f(xi*,ui*,t+i*dtau)*dtauを計算する
            prev_X_=x_;
            for i = 1:obj.N_step
                X_((obj.x_size*i)+1:obj.x_size*i+obj.x_size, 1)=prev_X_+obj.calModel(prev_X_, U_((obj.U_size*(i-1)+1):obj.U_size*(i-1)+obj.U_size, 1))*obj.dt;
                prev_X_=X_((obj.x_size*i)+1:obj.x_size*i+obj.x_size, 1);
            end
            obj.X_cal=X_;
            %4
            %lamda_(lamda*)を求める
            %lamdaN*=(rphi/rx)^T(xN*,t+T)を計算する
            Lamda_((obj.x_size*(obj.N_step-1))+1:(obj.x_size*(obj.N_step-1))+obj.x_size, 1)=(obj.rphirx(X_(obj.x_size*(obj.N_step-1)+1:obj.x_size*(obj.N_step-1)+obj.x_size, 1)))';
            %lamdai*=lamdai+1*+(rH/ru)^T*dtau
            prev_Lamda_=Lamda_(obj.x_size*(obj.N_step-1)+1:obj.x_size*(obj.N_step-1)+obj.x_size, 1);
            %逆順で解く
            %N_step-2の理由(N_step-1で最後のLamdaのグループなので(上でそこは計算してる),それの前だからN-2)
            for i = obj.N_step-2:-1:0
                Lamda_(obj.x_size*i+1:obj.x_size*i+obj.x_size, 1)=prev_Lamda_+(obj.rHrx(X_(obj.x_size*i+1:obj.x_size*i+obj.x_size, 1), U_(obj.U_size*i+1:obj.U_size*i+obj.U_size, 1), prev_Lamda_))'*obj.dt;
                prev_Lamda_=Lamda_(obj.x_size*i+1:obj.x_size*i+obj.x_size, 1);
            end
            %Fを求める
            for i = 0:obj.N_step-1
                lam_1=Lamda_((i*obj.x_size)+1, 1);
                lam_2=Lamda_((i*obj.x_size)+2, 1);
                lam_3=Lamda_((i*obj.x_size)+3, 1);
                u_1=U_((i*obj.U_size+1), 1);
                u_2=U_((i*obj.U_size+2), 1);
%                 dummy_=U_((i*obj.U_size+3), 1);
%                 mu_=U_((i*obj.U_size+4), 1);
                x_3=X_((i*obj.x_size+3), 1);
%                 F((i*obj.f_size)+1, 1)=obj.param_u_v*u_1+lam_1*cos(x_3)+lam_2*sin(x_3)+2*mu_*u_1;
%                 F((i*obj.f_size)+2, 1)=obj.param_u_omega*u_2+lam_3;
%                 F((i*obj.f_size)+3, 1)=obj.param_dummy*dummy_+2*mu_*obj.param_dummy*dummy_;
%                 F((i*obj.f_size)+4, 1)=u_1*u_1+dummy_*dummy_-100;
%                 F((i*obj.f_size)+1, 1)=obj.param_u_v*u_1+lam_1*cos(x_3)+lam_2*sin(x_3)+2*mu_*u_1;
                F((i*obj.f_size)+1, 1)=obj.param_u_v*u_1+lam_1*cos(x_3)+lam_2*sin(x_3);
                F((i*obj.f_size)+2, 1)=obj.param_u_omega*u_2+lam_3;
%                 F((i*obj.f_size)+3, 1)=obj.param_dummy*dummy_+2*mu_*obj.param_dummy*dummy_;
%                 F((i*obj.f_size)+4, 1)=u_1*u_1+dummy_*dummy_-100;
            end
        end
        function Av = calAv(obj, V_)
            U0= obj.U(1:obj.U_size,1);
            Av=(obj.calF(obj.U+(V_*obj.ht), obj.X+obj.calModel(obj.X, U0)*obj.ht)-obj.calF(obj.U, obj.X+obj.calModel(obj.X, U0)*obj.ht))/obj.ht;
        end
        function R0 = calR0(obj)
            %U'(0)=U0を使用する
            dX=obj.calModel(obj.X, obj.U(1:obj.u_size,1))*obj.ht;
%             disp("dX")
%             disp(dX)
%             disp("dX_end")
%             disp("U")
%             disp(obj.U)
%             disp("U_end")
%             disp("obj_F")
%             disp(obj.calF(obj.U, obj.X))
%             disp("obj_F_end")
%             disp("obj")
%             disp((obj.calF(obj.U, obj.X+dX)-obj.calF(obj.U, obj.X))/obj.ht-(obj.calF(obj.U+obj.U*obj.ht, obj.X+dX)-obj.calF(obj.U, obj.X+dX))/obj.ht)
%             disp("obj_end")
            R0=-1*obj.zeta*obj.calF(obj.U, obj.X)-(obj.calF(obj.U, obj.X+dX)-obj.calF(obj.U, obj.X))/obj.ht-(obj.calF(obj.U+obj.U*obj.ht, obj.X+dX)-obj.calF(obj.U, obj.X+dX))/obj.ht;
        end
    end
end      
