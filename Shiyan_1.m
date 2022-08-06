clear all
%实验数据 rms卡尔曼算法 加均方根噪声
rmsa=0.0848;
rmsd=0.00001;
ha=1/4096;
hd=1/4096;
nacce=load('C:\Users\Administrator\Desktop\\论文数据分析\论文钢梁数据\加速度计\acc1.txt');
ndisp=load('C:\Users\Administrator\Desktop\\论文数据分析\论文钢梁数据\视觉相机 高频\dis1.txt'); 
dt=ha; % acceleration sample interval
dd=hd;% displacemetn sample interval
sampletimes=hd/ha;% sample interval times.
mk=0;% times
x0=ndisp(1); % initial state estimate (displaement)
v0=0;% initial state estimate (velocity)

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

for t = 1:1:length(ndisp)*sampletimes-1
   
   u = nacce(t);
   Q_estimate = A * Q_estimate + B * u;
   Ppred = A * Pcorr * A' + Ex;
   loc_estimate_pred = [loc_estimate_pred; Q_estimate(1)];
   vel_estimate_pred = [vel_estimate_pred; Q_estimate(2)];
   Ppred_estimate = [Ppred_estimate;Ppred];
   Pcorr_estimate = [Pcorr_estimate;Ppred]; 
   loc_estimate = [loc_estimate; Q_estimate(1)];
%    vel_estimate = [vel_estimate; Q_estimate(2)]; 
   while t == round(mk*sampletimes+1)
    K = Ppred*C'* inv(C*Ppred*C'+ Ez./dd);
    Q_estimate = Q_estimate + K * (loc_meas(mk+1) - C * Q_estimate);
    Pcorr =  (eye(2)-K*C)*Ppred;
    mk = mk+1;
    Pcorr_estimate(2*t-1:2*t,:) = Pcorr;
    loc_estimate(t,1) = Q_estimate(1);
    vel_estimate(t,1) =  Q_estimate(2); 
   end  
    %% fixed-interval smoothing using RTS algorithm
    Ppred_estimate = [Ppred_estimate;Ppred];
    Pcorr_estimate = [Pcorr_estimate;Pcorr];
    Q_loc_estimate = [Q_loc_estimate; Q_estimate(1)];
    vel_estimate = [vel_estimate; Q_estimate(2)];
end

%%此后是作图程序，应该没什么问题
t_1(1)=0;
for n=1:1:length(Q_loc_estimate)-1
t_1(n+1)=t_1(n)+ha;
end
figure(1) %平滑后的位移图像
plot(t_1',Q_loc_estimate,'--');
hold on
t_2(1)=0;
for i=1:1:length(ndisp)-1
t_2(i+1)=t_2(i)+hd;
end
plot(t_2,ndisp);
% ts2(1)=0;
% for i=1:1:length(adis)-1
% ts2(i+1)=ts2(i)+ha;
% end
% plot(ts2,adis);
xlabel('time(sec)');ylabel('Disp');title('融合前后位移对比6');
legend('视觉位移','滤波结果')

% for j=1:length(disp)
%     rmss(j)=(dis(j)-Q_loc_estimate(j));
% end
% mean(abs(rmss))
% figure(2)
% t_3(1)=0;
% for i=1:1:length(Ndisp)-1
% t_3(i+1)=t_3(i)+hd;
% end
% plot(t_3,rmsD,'--');
% hold on
% plot(t_3,rmss');
% xlabel('time(sec)');ylabel('Disp error');title('加速度、位移噪声均为20%');
% legend('噪声','滤波误差')
rmss(1)=0;
for j=2:length(ndisp)-1
    rmss(j)=(ndisp(j)-Q_loc_estimate(j*sampletimes-4));
end
    rmsss=rms(rmss)
