clear all
%ʵ������ rms�������㷨 �Ӿ���������
rmsa=0.0848;
rmsd=0.00000001;
ha=1/4096;
hd=1/1024;
nacce=load('C:\Users\Administrator\Desktop\\�������ݷ���\���ĸ�������\���ٶȼ�\acc6.txt');
ndisp=load('C:\Users\Administrator\Desktop\\�������ݷ���\���ĸ�������\�Ӿ����������\dis6.txt'); 
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
    if N<=5 % N<5�������
       if t>N % ÿ��ǰһ�������һ���ͺ�ΪN�Ĺ⻬
           xsmooth = [loc_estimate(t-N:t,:) vel_estimate(t-N:t,:)]';% �⻬�����е�״̬��������
           xpred = [loc_estimate_pred(t-N:t,:)  vel_estimate_pred(t-N:t,:)]';% �⻬�����е�Э�������
           T=t; % �⻬��������Ϊ Pcorr_estimate��Ppred_estimate��ѭ������ 
           for k=N:-1:1 % Rnacce�⻬����
                F = Pcorr_estimate(2*T-3:2*T-2,:) * A' * inv(Ppred_estimate(2*T-1:2*T,:));   
                xsmooth(:,k) = xsmooth(:,k) + F * (xsmooth(:,k+1) - xpred(:,k+1));
                T=T-1; 
           end
           loc_smooth = [loc_smooth xsmooth(1,ceil(N/2))]; % �⻬���λ��
           vel_smooth = [vel_smooth xsmooth(2,ceil(N/2))]; % �⻬����ٶ� 
           if t==length(nacce) % �����ļ�������N��ֵ��һ������
              for k=N:-1:1
               loc_smooth = [loc_smooth loc_smooth(1,t-N)]; % �⻬���λ��
               vel_smooth = [vel_smooth vel_smooth(1,t-N)]; % �⻬����ٶ�  
              end
           end
       end
    elseif rem(length(nacce),N)==1 % ��N����nacce/N=1ʱ�Ĺ⻬����
        if length(nacce)-(m-1)*N~=1  % �ж�ʣ���ʱ���Ƿ���Դ���һ��N
          if t>m*N % ���в��ص���Rnacce�⻬
          xsmooth = [loc_estimate(t-N:t,:) vel_estimate(t-N:t,:)]'; % �⻬�����е�״̬��������
          xpred = [loc_estimate_pred(t-N:t,:)  vel_estimate_pred(t-N:t,:)]'; % �⻬�����е�Э�������
          T=t;
             for k=N:-1:1 % Rnacce�⻬����
               F = Pcorr_estimate(2*T-3:2*T-2,:) * A' * inv(Ppred_estimate(2*T-1:2*T,:));   
               xsmooth(:,k) = xsmooth(:,k) + F * (xsmooth(:,k+1) - xpred(:,k+1));
               T=T-1;
             end;
          loc_smooth = [loc_smooth xsmooth(1,1:N)]; % �⻬���λ��
          vel_smooth = [vel_smooth xsmooth(2,1:N)]; % �⻬����ٶ� 
          m=m+1;
          end   
        end
        if t==length(nacce) % �����һ��ʱ�̵�ֵ���н��ƹ���
          loc_smooth(1,t)= loc_smooth(1,t-1); % �⻬���λ��
          vel_smooth(1,t)= vel_smooth(1,t-1); % �⻬����ٶ� 
        end
    else % ��N>5��nacce/N������1�������
        if length(nacce)-(m-1)*N>N % �ж�ʣ���ʱ���Ƿ���Դ���һ��N��������Լ������У������ʣ�µĲ���һ��N��ʱ�̽���Rnacce�⻬
          if t>m*N % ���в��ص���Rnacce�⻬��ÿ��ǰN������һ�ι⻬
          xsmooth = [loc_estimate(t-N:t,:) vel_estimate(t-N:t,:)]';
          xpred = [loc_estimate_pred(t-N:t,:)  vel_estimate_pred(t-N:t,:)]';
          T=t; % �⻬��������Ϊ Pcorr_estimate��Ppred_estimate��ѭ������
          for k=N:-1:1 % Rnacce�⻬����
               F = Pcorr_estimate(2*T-3:2*T-2,:) * A' * inv(Ppred_estimate(2*T-1:2*T,:));   
               xsmooth(:,k) = xsmooth(:,k) + F * (xsmooth(:,k+1) - xpred(:,k+1));
               T=T-1;
          end; 
          loc_smooth = [loc_smooth xsmooth(1,1:N)]; % �⻬���λ��
          vel_smooth = [vel_smooth xsmooth(2,1:N)]; % �⻬����ٶ� 
          m=m+1;
          end  
        else
           if t==length(nacce) % ��ʣ�µĲ���һ��N��ʱ�̽���Rnacce�⻬
             xsmooth = [loc_estimate((m-1)*N+1:t,:) vel_estimate((m-1)*N+1:t,:)]'; % �⻬�����е�״̬��������
             xpred = [loc_estimate_pred((m-1)*N+1:t,:)  vel_estimate_pred((m-1)*N+1:t,:)]'; % �⻬�����е�Э�������
             T=length(nacce); % �⻬��������Ϊ Pcorr_estimate��Ppred_estimate��ѭ������
              for k=length(nacce)-(m-1)*N-1:-1:1
                F = Pcorr_estimate(2*T-3:2*T-2,:) * A' * inv(Ppred_estimate(2*T-1:2*T,:));   
                xsmooth(:,k) = xsmooth(:,k) + F * (xsmooth(:,k+1) - xpred(:,k+1));
                T=T-1; 
              end
             loc_smooth = [loc_smooth xsmooth(1,:)]; % �⻬���λ��
             vel_smooth = [vel_smooth xsmooth(2,:)]; % �⻬����ٶ� 
           end
        end
end
end
% output displacment and velocity estimate
nacce=nacce';ndisp=ndisp';loc_smooth=loc_smooth';vel_smooth=vel_smooth';
%% figure1 �˲���λ����ԭʼλ�ƱȽ�
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
xlabel('time(sec)');ylabel('Disp');title('���ٶȡ�λ��������Ϊ10%');
legend('��ʵλ��','�˲����')



% %% figure2 �˲�����������ľ����������ıȽ�
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
% xlabel('time(sec)');ylabel('Disp error');title('���ٶȡ�λ��������Ϊ10%');
% legend('����','�˲����');
%% figure3 ƽ����λ�����˲�λ�ƱȽ�
figure(3)
plot(t_1',loc_estimate);
hold on
t_5(1)=0;
for i=1:1:length(loc_smooth)-1
t_5(i+1)=t_5(i)+ha;
end
plot(t_5,loc_smooth','--');
xlabel('time(sec)');ylabel('Disp');title('���ٶȡ�λ��������Ϊ10%');
legend('�˲�λ��','ƽ��λ��');
%% figure4 ƽ����λ��������˲�λ�����ıȽ�
% figure(4)
% for j=1:length(loc_smooth)
%     rms_smooth(j)=(dis(j)-loc_smooth(j));
% end
% rms_loc_smooth=sum(rms_smooth.*rms_smooth);
% rms_loc_smooth=(rms_loc_smooth/length(rms_loc_smooth))^0.5
% plot(t_1,rms_estimate');
% hold on
% plot(t_1,rms_smooth','--');
% xlabel('time(sec)');ylabel('Disp error');title('���ٶȡ�λ��������Ϊ10%');
% legend('�˲����','ƽ�����');
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