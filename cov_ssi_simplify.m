function [freq,modal_shape,damping] = cov_ssi_simplify()
% cov_ssi����ʵ��ģ̬ʶ��ĳ��򣨻��������µ�ģ̬ʶ��
% ���������ڳ����ڲ����У���Ҫͨ���������ݱ���ʵ��
% �����freq,modal_shape,damping ����ΪƵ�ʡ����ͺ�����Ⱦ���������
% ��������
% ���������д����ڶ���excel����ʱ����Ҫ���ݸ�ʽΪ����һ��Ϊʱ�����У���λΪ�룬֮��ÿһ�б�ʾһ�����ļ��ٶȡ��ٶȻ���λ������
[FileName,PathName]=uigetfile('*.xlsx','��ѡ��Ҫ�������ݵ��Ӷ�ʱ�̱��');  %��ʾѡ��ʱ���ļ��ĶԻ���
File=strcat(PathName,FileName);  %�����������ӳɵ����ַ���
time_history_0 = xlsread(File);   % �ò�����Ҫ���ݸ�ʽΪ����һ��Ϊʱ�����У���λΪ�룬֮��ÿһ�б�ʾһ�����ļ��ٶȡ��ٶȻ���λ������
[m,n] = size(time_history_0);
if m < 10 || n < 2
    key_of_source = questdlg('��⵽�������ݳ��ȹ�С��û�в�����ݣ��Ƿ����ɽ��в���ʶ��','�Ƿ����ɽ��в���ʶ��','��','��','��');
    if key_of_source == '��'
        time_history = time_history_0;
    else
        freq = [];modal_shape = [];damping = [];
        return
    end
else
    time_history = time_history_0;
end

T = time_history(3,1) - time_history(2,1);    % �������ʱ����
fs = 1/T;         % �������Ƶ��

% cov_ssi ģ̬ʶ�����������Դ�����磬ʱ��̫�������ǳ����ˣ��������
[TimePointNum,MesNodeTotalNum]=size(time_history_0);
MesNodeTotalNum=MesNodeTotalNum-1;

prompt={'TOEPLITZ�������M(>=����Ƶ��/(2*�ṹ��Ƶ)):','ϵͳ����N(<i*ʵ�ʲ����):'};
dlg_title='������COV_SSI��������ı�Ҫ����(���Ϊ2��������)';
num_lines=[1 100;1 100];
def={'50','30'};                                 %ϵͳ����N=2*n
answer=inputdlg(prompt,dlg_title,num_lines,def,'on');

M=str2double(answer{1});
N=str2double(answer{2});
disp('----------ģ̬����ʶ��ʼ----------');
disp('          ����Hankel������...');
pause(0.000001);  %�ó�����ͣһ���ʱ�䣬��ˢ�½���
j=TimePointNum-2*M;
Hankel=zeros((2*M+1)*MesNodeTotalNum,j);
for ti=1:j  %���λ����forѭ����Ч�ʱȽϵͣ���������ʵ�־���Ĺ�������Ч�ʻ��һ��
    m=0;
    for tii=1:(2*M+1)
        km=ti+tii-1;
        m=m+1;
        %����ÿ�θ�ֵ�������е���ʼλ��
        beginy=(m-1)*MesNodeTotalNum+1;
        %����ÿ�θ�ֵ�������еĽ���λ��
        endy=m*MesNodeTotalNum;
        Hankel(beginy:endy,ti)=(time_history_0(km,2:(MesNodeTotalNum+1)))';
    end
end

Yp=Hankel(1:M*MesNodeTotalNum,:);
Yf=Hankel(M*MesNodeTotalNum+1:2*M*MesNodeTotalNum,:);
Yf2=Hankel((M+1)*MesNodeTotalNum+1:(2*M+1)*MesNodeTotalNum,:);
%��ɵ�һ��Toeplitz����
Teop1=Yf*Yp'/j;
%��ɵڶ���Toeplitz����
Teop2=Yf2*Yp'/j;
disp('          ���ڽ�������ֵ�ֽ�...');
pause(0.000001);  %�ó�����ͣһ���ʱ�䣬��ˢ�½���
[U,S,V]=svd(Teop1);

%��ȡǰN������ֵ�ֽ���
U_de=U(:,1:N);
S_de=S(1:N,1:N);
V_de=V(:,1:N);

%����C,A����Ĺ���
O_GJ=U_de*S_de^0.5;
C_GJ=O_GJ(1:MesNodeTotalNum,:);
A_GJ=S_de^-0.5*(U_de')*Teop2*V_de*S_de^-0.5;

disp('         ����ֵ�ֽ���ɣ�����ģ̬����ʶ��...');
pause(0.000001);  %�ó�����ͣһ���ʱ�䣬��ˢ�½���
%����״̬����A������ֵ����������
[V_of_A,D_of_A]=eig(A_GJ);


%����ϵͳ������ֵ����������
% fs=user.data.sim_freq;
% T=1/user.data.sim_freq;
Total_mode_jieShu=N;
%��������
MODE=C_GJ*V_of_A;
%����Ƶ��
LN_lamuda = zeros(N,1);
Re_Abs = zeros(N,1);
Im_Abs = zeros(N,1);
Freq = zeros(N,1);
Dpro = zeros(N,1);
for i = 1:N
    LN_lamuda(i,1)=log(D_of_A(i,i))/T;
    Re_Abs(i,1)=abs(real(LN_lamuda(i,1)));
    Im_Abs(i,1)=abs(imag(LN_lamuda(i,1)));
    Freq(i,1)=abs(LN_lamuda(i,1))/2/pi;      %�������õ���Ƶ�ʵ�λΪHZ
    %���������
    Dpro(i,1)=Re_Abs(i,1)/abs(LN_lamuda(i,1));
end
disp('          ʶ�����б���...');
pause(0.000001);  %�ó�����ͣһ���ʱ�䣬��ˢ�½���
%��ʵ�����鲽�ֱ��һ��
for opp_index=1:Total_mode_jieShu
    min_of_real=min(real(MODE(:,opp_index)));
    max_of_real=max(real(MODE(:,opp_index)));
    min_of_imag=min(imag(MODE(:,opp_index)));
    max_of_imag=max(imag(MODE(:,opp_index)));
    if max_of_real~=min_of_real
        Mode_value_Total_of_real_GYH=2*(real(MODE(:,opp_index))-min_of_real)/(max_of_real-min_of_real)-1;
    else
        Mode_value_Total_of_real_GYH=zeros(MesNodeTotalNum,1);
    end
    if max_of_imag~=min_of_imag
        Mode_value_Total_of_imag_GYH=2*(imag(MODE(:,opp_index))-min_of_imag)/(max_of_imag-min_of_imag)-1;
    else
        Mode_value_Total_of_imag_GYH=zeros(MesNodeTotalNum,1);
    end
    MODE(:,opp_index)=Mode_value_Total_of_real_GYH+Mode_value_Total_of_imag_GYH*sqrt(-1);
    if imag(MODE(:,opp_index))~=0
        MODE(:,opp_index)=imag(MODE(:,opp_index));
    else
        MODE(:,opp_index)=real(MODE(:,opp_index));
    end
end

% ��ʵģ̬�ж�
% ��������޳����ģ̬�����ǲ�֪������ȴ�Χ����������һ���Ƚϴ��ֵ
prompt={'��������ʵģ̬��������ֵ�������޳������ģ̬����:'};  %j>20M
dlg_title='Input Parameter of Calculate';
num_lines=1;
def={'0.05'};
answer=inputdlg(prompt,dlg_title,num_lines,def);
user.max_of_damping=str2double(answer{1,1});
M2=length(Freq);
mm2=1;
for m2=1:M2
    if (Freq(m2)<=(0.5*fs))&&(Freq(m2)>0)&&(Dpro(m2)>0)&&(Dpro(m2)<user.max_of_damping)      % ����0.05��Ҫ���û�����
        Freq2(mm2)=Freq(m2);
        Dpro2(mm2)=Dpro(m2);
        MODE2(:,mm2)=MODE(:,m2);
        mm2=mm2+1;
    end
end

[Freq3,IX]=sort(Freq2);
for m3=1:length(Freq3)
    Dpro3(m3)=Dpro2(IX(m3));
    MODE3(:,m3)=MODE2(:,IX(m3));
end

% ��Ƶ�ж�ϵ��
crtifreq = 1;              % Ƶ���жϱ�׼
crtidamping = 5;           % ������жϱ�׼
crtimodeshape = 99.9;      % ģ̬�жϱ�׼��MAC��
M3=length(Freq3);
dualfreq=ones(1,M3);
for m3=1:M3
    for n3=m3+1:M3
        X=MODE3(:,m3);
        Y=MODE3(:,n3);
        MAC_shape=(X'*Y)^2/((X'*X)*(Y'*Y));             % ģ̬��֤׼�򣬱�������ģ̬�������
        if (((Freq3(n3)-Freq3(m3))*100/Freq3(m3))<crtifreq)...
                && (abs(((Dpro3(m3)-Dpro3(n3))*100/Dpro3(m3))) < crtidamping)...
                && (abs(MAC_shape*100) > crtimodeshape)
            dualfreq(n3)=0;
        end
    end
end
indicefreq=find(dualfreq);
M4=length(indicefreq);
for m4=1:M4
    Freq4(m4)=Freq3(indicefreq(m4));
    Dpro4(m4)=Dpro3(indicefreq(m4));
    MODE4(:,m4)=MODE3(:,indicefreq(m4));
end
disp('----------ʶ�����----------');
pause(0.000001);  %�ó�����ͣһ���ʱ�䣬��ˢ�½���


% ģ̬����ʶ�����б���ʾ
user.results.freq=Freq4; user.results.damping=Dpro4; user.results.mode_shape=MODE4;
% ���ȫ������Ƶ��
[~,n]=size(user.results.freq);
c_n={('����'),('Ƶ��(Hz)'),('�����(%)')};
b=zeros(n,3);
for i=1:n
    b(i,1)=i;
    b(i,2)=Freq4(i);
    b(i,3)=100*Dpro4(i);
end
uitable('Tag', 'out_table',...
    'Units', 'Normalized',...
    'Data',b,...
    'ColumnName',c_n,...
    'FontSize',12,...
    'Position', [0.02,0.15,0.8,0.8]);
pause(0.000001);  %�ó�����ͣһ���ʱ�䣬��ˢ�½���
freq = user.results.freq;
damping = user.results.damping;
modal_shape = user.results.mode_shape;

% ����Ϊ�ȶ�ͼ��������Ƴ���
key_of_source = questdlg('��⵽ʶ���Ѿ���ɣ��Ƿ�����ȶ�ͼ���ƣ�','�Ƿ�����ȶ�ͼ���ƣ�','��','��','��');
if key_of_source == '��'
    return
else
    p = 1;        % �ȶ�ͼ���ݴ洢�������м�������
    N_min = 2;     % �����ȶ�ͼ����С����
    if rem(N,2) == 0
        N_max = N;    % �����ȶ�ͼ��������
    else
        N_max = N+1;    % �����ȶ�ͼ��������
    end
    
    S_freq=-ones((N_max-N_min)/2+1,N_max);
    S_damp=-ones((N_max-N_min)/2+1,N_max);
    
    for N_system=N_min:2:N_max
        U_N=U(:,1:N_system);
        V_N=V(:,1:N_system);
        S_N=S(1:N_system,1:N_system);
        A_N=S_N^-0.5*(U_N')*Teop2*V_N*S_N^-0.5;
        [~,D_N]=eig(A_N);
        [m,~]=size(D_N);
        Eignvalue = zeros(1,m);
        for i=1:m
            Eignvalue(1,i)=D_N(i,i);
        end
        Eignvalue_Ac=log(Eignvalue)/T;
        Freq_N=abs(Eignvalue_Ac)/2/pi;
        Damp_N=-real(Eignvalue_Ac)./abs(Eignvalue_Ac);
        for q=1:N_system
            S_freq(p,q)=Freq_N(1,q);
            S_damp(p,q)=Damp_N(1,q);
        end
        p=p+1;
    end
    
    
    % �������֮ǰ��ע�͵��ģ�2021��4��28���޸ĺ�ȡ��ע��,������ȿ������޳����ģ̬
    [b_max,a_max] = size(S_damp);
    try
        for a=1:a_max
            for b=1:b_max
                if S_damp(b,a)>=user.max_of_damping || S_damp(b,a)<=0
                    S_freq(b,a)=-100;
                end
            end
        end
    catch
        prompt={'��������ʵģ̬��������ֵ�������޳������ģ̬����:'};  %j>20M
        dlg_title='Input Parameter of Calculate';
        num_lines=1;
        def={'0.05'};
        answer=inputdlg(prompt,dlg_title,num_lines,def);
        user.max_of_damping=str2double(answer{1,1});
        for a=1:a_max
            for b=1:b_max
                if S_damp(b,a)>=user.max_of_damping || S_damp(b,a)<=0
                    S_freq(b,a)=-100;
                end
            end
        end
    end
    
    
    % *******************************�����ȶ�ͼ********************************
    [m,n] = size(time_history_0);
    y_axis=N_min:2:N_max;        % �ȶ�ͼy�������
    while true
        prompt={'������PSD�׼���ʱ�Ĵ�����:'};
        dlg_title='Input Parameter of Calculate';
        num_lines=1;
        def={'1024'};
        answer=inputdlg(prompt,dlg_title,num_lines,def);
        correction=str2double(answer{1,1});
        if correction > 0 && 2*correction < m
            nfft = correction;
            break
        else
            disp('          PSD�׼���ʱ�Ĵ���������ֵ����');
            disp('          ����������С��ʶ�����ݳ���һ�����������');
        end
    end
    
    [Pxx,f]=pwelch(time_history_0(:,2:n),hanning(nfft),nfft/2,nfft,fs);
    pxx = mean(Pxx,2);
    figure0 = figure;
    hp = axes('parent',figure0,...
        'FontSize',10,...
        'Units', 'Normalized',...
        'Nextplot', 'Add',...
        'Box', 'on',...
        'Position',[0.15 0.15 0.8 0.8]);
    yyaxis(hp,'left')
    plot(hp,S_freq,y_axis,'*');
    set(hp,'xlim',[0 fs/2]);   %����������̶�
    set(hp,'ylim',[0 N_max]);
    set(hp.XLabel,'string','Frequency(Hz)','FontSize',10);   %�����������ǩ
    set(hp.YLabel,'string','Number of Singular Values(N)','FontSize',10);
    yyaxis(hp,'right')
    plot(hp,f,pxx,'-');
    set(hp,'xlim',[0 fs/2]);   %����������̶�
    set(hp.XLabel,'string','Frequency(Hz)','FontSize',10);   %�����������ǩ
    set(hp.YLabel,'string','Power Spectral Density(dB)','FontSize',10);
    grid on
    title('�ȶ�ͼ','FontSize',10);
    disp('----------�ȶ�ͼ������ɣ�----------')
end

end