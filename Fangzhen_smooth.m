clear all
%�������� �ɵĿ������㷨 �Ӿ��������� RTS�ͺ�5ƽ��
N=100;
a=-0.02;b=0.1;c=0.5;%�����ʼ����
n=1;acce(1)=2;nacce2(1)=0;vol(1)=0.9;disp(1)=0;nacce1(1)=0;dis(1)=0;nacce3(1)=0;
vol(2)=1+(-0.1+2*0.002)*cos((-0.1+0.002)*0.002);
while nacce2<20-0.001  %���ٶȼ��㣬����Ƶ��Ϊ1000
nacce2(n+1)=nacce2(n)+0.001;
n=n+1;
acce(n)=2*c*cos((a+c*nacce2(n))*nacce2(n))-(a+2*c*nacce2(n))^2*sin((a+c*nacce2(n))*nacce2(n));
dis(n)=sin((a+c*nacce2(n))*nacce2(n))+b*nacce2(n);
vol(n)=b+(a+2*c*nacce2(n))*cos((a+c*nacce2(n))*nacce2(n));
end
n=1;
while nacce1<20-0.01  %λ�Ƽ��㣬����Ƶ��Ϊ100
nacce1(n+1)=nacce1(n)+0.01;
n=n+1;
disp(n)=sin((a+c*nacce1(n))*nacce1(n))+b*nacce1(n);
end

ha=1/1000;
hd=1/100;

% Nacce=load('acc.txt');
% Ndisp=load('dis1.txt');
% dis=load('dis.txt');
% adis=dis(100:149,1);
% acce=Nacce(496:746,1);
% disp=Ndisp(100:149,1);
% 
% rmsA=randn(1,length(acce)); 
% rmsA=rmsA/std(rmsA); 
% rmsA=rmsA-mean(rmsA); 
% rmsa=sqrt(0.1); 
% rmsA=rmsA*rmsa;
% nacce=acce+rmsA';
% 
% rmsD=randn(1,length(disp)); 
% rmsD=rmsD/std(rmsD); 
% rmsD=rmsD-mean(rmsD); 
% rmsd=sqrt(0.1); 
% rmsD=rmsD*rmsd;
% ndisp=disp+rmsD';
% %rmsd=0.0000001;

rmsa=sum(acce.*acce);
rmsa=(rmsa/length(acce))^0.5;
nacce=acce+0.1*rmsa*randn(1,length(acce));%   ���ٶ�����

rmsd=sum(disp.*disp);
rmsd=(rmsd/length(disp))^0.5;
ndisp=disp+0.1*rmsd*randn(1,length(disp));%   λ������
%% ����ȡ��ʱ�� Ƶ��
dt = ha;             % ���ٶȲ���ʱ��
dd = hd;             % λ�Ʋ���ʱ��
sampletimes = hd/ha; % �������
mk = 0.0;            % ȡ���ݴ���
x0 = 0.0;            % ��ʼλ��
v0 = 0.0;            % ��ʼ�ٶ�


%% ������¾���
A = [1 dt; 0 1] ; % state transition matrix:  expected flight of the Quail (state prediction)
B = [dt^2/2; dt]; % input control matrix:  expected effect of the input accceleration on the state.
C = [1 0];        % measurement matrix: the expected measurement given the predicted state (likelihood)
%since we are only measuring position (too hard for the ninja to calculate speed), we set the velocity variable to
%zero.

%% ����������
Q= [x0; v0];      % initized state--it has two componennacce: [position; velocity] 
Q_estimate = Q;   % x_estimate of initial location estimation of where the Quail is (what we are updating)
QuailAccel_noise_mag =var(nacce-acce);% process noise
NinjaVision_noise_mag =var(ndisp-disp); % measurement noise
Ex = QuailAccel_noise_mag* [dt^3/3 dt^2/2; dt^2/2 dt]; % Ex convert the process noise (stdv) into covariance matrix
Ez = NinjaVision_noise_mag ;  % Ez convert the measurement noise (stdv) into covariance matrix
Pcorr = eye(2);%Ex;  % ;  %estimate of initial Quail position variance (covariance matrix)
loc_smooth = [];% smooth of displacment
vel_smooth = [];% smooth of velocity

%% ��ʼ���û��������

loc_meas = ndisp; 

%% ������������
%��ʼ�����Ʊ���
loc_estimate = []; % λ�ƹ���
vel_estimate = [];   % �ٶȹ���
loc_estimate_pred = [];
vel_estimate_pred = [];
Pcorr_estimate = [];
Ppred_estimate = [];
Q = [x0; v0]; % ���³�ʼ��״̬
m=1; % N>5������½�����ʱ��ֳ�m�ηֱ����Rnacce�⻬

for t = 1:1:length(nacce)
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
disp=disp';dis=dis';nacce=nacce';ndisp=ndisp';loc_smooth=loc_smooth';vel_smooth=vel_smooth';acce=acce';

%% figure1 �˲���λ����ԭʼλ�ƱȽ�
t_1(1)=0;
for n=1:1:length(loc_estimate)-1
t_1(n+1)=t_1(n)+ha;
end
figure(1) 
plot(t_1',loc_estimate,'--');
hold on
t_2(1)=0;
for i=1:1:length(disp)-1
t_2(i+1)=t_2(i)+hd;
end
plot(t_2,disp);
xlabel('time(sec)');ylabel('Disp');title('���ٶȡ�λ��������Ϊ10%');
legend('��ʵλ��','�˲����')



%% figure2 �˲�����������ľ����������ıȽ�
for j=1:length(dis)
    rms_estimate(j)=(dis(j)-loc_estimate(j));
end
rms_loc_extimate=sum(rms_estimate.*rms_estimate);
rms_loc_extimate=(rms_loc_extimate/length(rms_loc_extimate))^0.5
figure(2)
t_3(1)=0;
for i=1:1:length(ndisp)-1
t_3(i+1)=t_3(i)+hd;
end
plot(t_3,ndisp-disp,'--');
hold on
plot(t_1,rms_estimate','Color','black','LineWidth',2);
xlabel('time(sec)');ylabel('Disp error');title('���ٶȡ�λ��������Ϊ10%');
legend('����','�˲����');
%% figure3 ƽ����λ�����˲�λ�ƱȽ�
figure(3)
plot(t_1',loc_estimate);
hold on
plot(t_1,loc_smooth','--');
xlabel('time(sec)');ylabel('Disp');title('���ٶȡ�λ��������Ϊ10%');
legend('�˲�λ��','ƽ��λ��');
%% figure4 ƽ����λ��������˲�λ�����ıȽ�
figure(4)
for j=1:length(loc_smooth)
    rms_smooth(j)=(dis(j)-loc_smooth(j));
end
rms_loc_smooth=sum(rms_smooth.*rms_smooth);
rms_loc_smooth=(rms_loc_smooth/length(rms_loc_smooth))^0.5
plot(t_1,rms_estimate');
hold on
plot(t_1,rms_smooth','--');
xlabel('time(sec)');ylabel('Disp error');title('���ٶȡ�λ��������Ϊ10%');
legend('�˲����','ƽ�����');

% figure(5)
% t_4(1)=0;
% for n=1:1:length(loc_estimate)-1
%     t_4(n+1)=t_4(n)+1;
% end
% plot(t_4,loc_estimate','--');
% hold on
% plot(t_4,dis');

