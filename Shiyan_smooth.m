clear all
%实验数据 rms卡尔曼算法 加均方根噪声
rmsa=0.0848;
rmsd=0.00000001;
ha=1/4096;
hd=1/1024;
nacce=load('C:\Users\Administrator\Desktop\\论文数据分析\论文钢梁数据\加速度计\acc6.txt');
ndisp=load('C:\Users\Administrator\Desktop\\论文数据分析\论文钢梁数据\视觉相机降采样\dis6.txt'); 
N=5;

dt=ha; 
dd=hd;
sampletimes=hd/ha;
mk=0;
x0=ndisp(1); 
v0=0;

A = [1 dt; 0 1] ; 
B = [dt^2/2; dt]; 
C = [1 0];       

Q_loc_estimate = [];% smooth of displacment
vel_smooth = [];% smooth of velocity
Q_estimate=[x0;v0];% x_estimate of initial location estimation
QuailAccel_noise_mag =rmsa;% process noise
NinjaVision_noise_mag =rmsd; % measurement noise
Ex = QuailAccel_noise_mag^2 * [dt^3/3 dt^2/2; dt^2/2 dt];
Ez = NinjaVision_noise_mag^2; 
Pcorr = [1 0;1 0];

loc_meas = ndisp;% displacment measurement

loc_estimate = []; % displacment estimate
vel_estimate = [];  % velocity estimate  
loc_estimate_pred = [];
vel_estimate_pred = [];
Pcorr_estimate = [];
Ppred_estimate = [];
Q = [x0; v0]; 
loc_smooth = [];% smooth of displacment
vel_smooth = [];% smooth of velocity
m=1;

for t = 1:1:length(ndisp)*sampletimes-1
    %predict next state of the quail with the last state and predicted motion.
    u = nacce(t);
    Q_estimate = A * Q_estimate + B * u;   
    %predict next covariance
    Ppred = A * Pcorr * A' + Ex;        
    loc_estimate_pred = [loc_estimate_pred; Q_estimate(1)];
    vel_estimate_pred = [vel_estimate_pred; Q_estimate(2)];
    
    Ppred_estimate = [Ppred_estimate;Ppred];
    Pcorr_estimate = [Pcorr_estimate;Ppred]; 
    loc_estimate = [loc_estimate; Q_estimate(1)];
    vel_estimate = [vel_estimate; Q_estimate(2)]; 
    
    %multirate data fusion
    while t == round(mk*sampletimes+1)
    % Kalman Gain
    K = Ppred*C'* inv(C*Ppred*C'+ Ez./dd);
    % Update the state estimate.
    Q_estimate = Q_estimate + K * (loc_meas(mk+1) - C * Q_estimate);
    % update covariance estimation.
    Pcorr =  (eye(2)-K*C)*Ppred;
    mk = mk+1;
    %Store for plotting
    Pcorr_estimate(2*t-1:2*t,:) = Pcorr;
    loc_estimate(t,1) = Q_estimate(1);
    vel_estimate(t,1) =  Q_estimate(2);
    end
   %% fixed-interval smoothing using Rnacce algorithm
    if N<=5 % N<5的情况下
       if t>N % 每向前一步便进行一次滞后为N的光滑
           xsmooth = [loc_estimate(t-N:t,:) vel_estimate(t-N:t,:)]';% 光滑过程中的状态向量矩阵
           xpred = [loc_estimate_pred(t-N:t,:)  vel_estimate_pred(t-N:t,:)]';% 光滑过程中的协方差矩阵
           T=t; % 光滑过程中作为 Pcorr_estimate和Ppred_estimate的循环变量 
           for k=N:-1:1 % Rnacce光滑过程
                F = Pcorr_estimate(2*T-3:2*T-2,:) * A' * inv(Ppred_estimate(2*T-1:2*T,:));   
                xsmooth(:,k) = xsmooth(:,k) + F * (xsmooth(:,k+1) - xpred(:,k+1));
                T=T-1; 
           end
           loc_smooth = [loc_smooth xsmooth(1,ceil(N/2))]; % 光滑后的位移
           vel_smooth = [vel_smooth xsmooth(2,ceil(N/2))]; % 光滑后的速度 
           if t==length(nacce) % 把最后的几个不足N的值做一个近似
              for k=N:-1:1
               loc_smooth = [loc_smooth loc_smooth(1,t-N)]; % 光滑后的位移
               vel_smooth = [vel_smooth vel_smooth(1,t-N)]; % 光滑后的速度  
              end
           end
       end
    elseif rem(length(nacce),N)==1 % 当N满足nacce/N=1时的光滑过程
        if length(nacce)-(m-1)*N~=1  % 判断剩余的时刻是否可以凑足一个N
          if t>m*N % 进行不重叠的Rnacce光滑
          xsmooth = [loc_estimate(t-N:t,:) vel_estimate(t-N:t,:)]'; % 光滑过程中的状态向量矩阵
          xpred = [loc_estimate_pred(t-N:t,:)  vel_estimate_pred(t-N:t,:)]'; % 光滑过程中的协方差矩阵
          T=t;
             for k=N:-1:1 % Rnacce光滑过程
               F = Pcorr_estimate(2*T-3:2*T-2,:) * A' * inv(Ppred_estimate(2*T-1:2*T,:));   
               xsmooth(:,k) = xsmooth(:,k) + F * (xsmooth(:,k+1) - xpred(:,k+1));
               T=T-1;
             end;
          loc_smooth = [loc_smooth xsmooth(1,1:N)]; % 光滑后的位移
          vel_smooth = [vel_smooth xsmooth(2,1:N)]; % 光滑后的速度 
          m=m+1;
          end   
        end
        if t==length(nacce) % 对最后一个时刻的值进行近似估计
          loc_smooth(1,t)= loc_smooth(1,t-1); % 光滑后的位移
          vel_smooth(1,t)= vel_smooth(1,t-1); % 光滑后的速度 
        end
    else % 在N>5且nacce/N不等于1的情况下
        if length(nacce)-(m-1)*N>N % 判断剩余的时刻是否可以凑足一个N，如果可以继续运行，否则对剩下的不足一个N的时刻进行Rnacce光滑
          if t>m*N % 进行不重叠的Rnacce光滑，每向前N步进行一次光滑
          xsmooth = [loc_estimate(t-N:t,:) vel_estimate(t-N:t,:)]';
          xpred = [loc_estimate_pred(t-N:t,:)  vel_estimate_pred(t-N:t,:)]';
          T=t; % 光滑过程中作为 Pcorr_estimate和Ppred_estimate的循环变量
          for k=N:-1:1 % Rnacce光滑过程
               F = Pcorr_estimate(2*T-3:2*T-2,:) * A' * inv(Ppred_estimate(2*T-1:2*T,:));   
               xsmooth(:,k) = xsmooth(:,k) + F * (xsmooth(:,k+1) - xpred(:,k+1));
               T=T-1;
          end; 
          loc_smooth = [loc_smooth xsmooth(1,1:N)]; % 光滑后的位移
          vel_smooth = [vel_smooth xsmooth(2,1:N)]; % 光滑后的速度 
          m=m+1;
          end  
        else
           if t==length(nacce) % 对剩下的不足一个N的时刻进行Rnacce光滑
             xsmooth = [loc_estimate((m-1)*N+1:t,:) vel_estimate((m-1)*N+1:t,:)]'; % 光滑过程中的状态向量矩阵
             xpred = [loc_estimate_pred((m-1)*N+1:t,:)  vel_estimate_pred((m-1)*N+1:t,:)]'; % 光滑过程中的协方差矩阵
             T=length(nacce); % 光滑过程中作为 Pcorr_estimate和Ppred_estimate的循环变量
              for k=length(nacce)-(m-1)*N-1:-1:1
                F = Pcorr_estimate(2*T-3:2*T-2,:) * A' * inv(Ppred_estimate(2*T-1:2*T,:));   
                xsmooth(:,k) = xsmooth(:,k) + F * (xsmooth(:,k+1) - xpred(:,k+1));
                T=T-1; 
              end
             loc_smooth = [loc_smooth xsmooth(1,:)]; % 光滑后的位移
             vel_smooth = [vel_smooth xsmooth(2,:)]; % 光滑后的速度 
           end
        end
end
end
% output displacment and velocity estimate
nacce=nacce';ndisp=ndisp';loc_smooth=loc_smooth';vel_smooth=vel_smooth';
%% figure1 滤波的位移与原始位移比较
t_1(1)=0;
for n=1:1:length(loc_estimate)-1
t_1(n+1)=t_1(n)+ha;
end
figure(1) 
plot(t_1',loc_estimate,'--');
hold on
t_2(1)=0;
for i=1:1:length(ndisp)-1
t_2(i+1)=t_2(i)+hd;
end
plot(t_2,ndisp);
xlabel('time(sec)');ylabel('Disp');title('加速度、位移噪声均为10%');
legend('真实位移','滤波结果')



% %% figure2 滤波的误差与加入的均方根噪声的比较
% for j=1:length(dis)
%     rms_estimate(j)=(dis(j)-loc_estimate(j));
% end
% rms_loc_extimate=sum(rms_estimate.*rms_estimate);
% rms_loc_extimate=(rms_loc_extimate/length(rms_loc_extimate))^0.5
% figure(2)
% t_3(1)=0;
% for i=1:1:length(ndisp)-1
% t_3(i+1)=t_3(i)+hd;
% end
% plot(t_3,ndisp-disp,'--');
% hold on
% plot(t_1,rms_estimate','Color','black','LineWidth',2);
% xlabel('time(sec)');ylabel('Disp error');title('加速度、位移噪声均为10%');
% legend('噪声','滤波误差');
%% figure3 平滑的位移与滤波位移比较
figure(3)
plot(t_1',loc_estimate);
hold on
t_5(1)=0;
for i=1:1:length(loc_smooth)-1
t_5(i+1)=t_5(i)+ha;
end
plot(t_5,loc_smooth','--');
xlabel('time(sec)');ylabel('Disp');title('加速度、位移噪声均为10%');
legend('滤波位移','平滑位移');
%% figure4 平滑的位移误差与滤波位移误差的比较
% figure(4)
% for j=1:length(loc_smooth)
%     rms_smooth(j)=(dis(j)-loc_smooth(j));
% end
% rms_loc_smooth=sum(rms_smooth.*rms_smooth);
% rms_loc_smooth=(rms_loc_smooth/length(rms_loc_smooth))^0.5
% plot(t_1,rms_estimate');
% hold on
% plot(t_1,rms_smooth','--');
% xlabel('time(sec)');ylabel('Disp error');title('加速度、位移噪声均为10%');
% legend('滤波误差','平滑误差');
% 
% % figure(5)
% % t_4(1)=0;
% % for n=1:1:length(loc_estimate)-1
% %     t_4(n+1)=t_4(n)+1;
% % end
% % plot(t_4,loc_estimate','--');
% % hold on
% % plot(t_4,dis');
% 
rms_estimate(1)=0;
for j=2:length(ndisp)-1
    rms_estimate(j)=(ndisp(j)-loc_estimate(j*sampletimes-4));
end
    rmss_estimate=rms(rms_estimate)
    
rms_smooth(1)=0;
for j=2:length(ndisp)-1
    rms_smooth(j)=(ndisp(j)-loc_smooth(j*sampletimes-6));
end
    rmss_smooth=rms(rms_smooth)
 [i1,j1]=max(rms_estimate(1,140:end));
 wucha=i1/ndisp(j1)