clc;clear;
close all;
clear funP
import Config.*
global Xresult; % SSA
Xresult = [];
global Xresult2;    % PSO
Xresult2 = [];
global Xresult3;    % SSA_opt
Xresult3 = [];

load w

%% ���߻�������
f = 3.5e9; % Ƶ��
N = Config.N; % ��Ԫ����
dx = 3e8/f/2; % ��Ԫ���Ϊ������һ��
dy = 0.043;  %��Ϊy����ֻ��1�����ߵ�Ԫ����ֵ����ν
%% ��������λ��
Amp = ones(N,1);
Amp = abs(w*Amp);
% Phase = [0:N-1]'*10;
% Phase = ones(N,1) * 90;
% Phase = rand(N,1) * 180; % ����
Phase =  randn(N,1)*10 + 90; % ��̬

[F_dB] = PlanerArrayPattern(f,Amp,Phase,dx,dy);
% [HPBW,SLL] = Get_HPBW_and_SLL(0,F_dB);
[HPBW,SLL,Nulls,~] = Get_HPBW_SLL_Nulls(0,F_dB);
%% ���õ�ԪʧЧ
offarray=[18 19];
Amp2 = Amp;
Amp2(offarray,1) = 0;
[F_dB2] = PlanerArrayPattern(f,Amp2,Phase,dx,dy);
% [HPBW2,SLL2] = Get_HPBW_and_SLL(0,F_dB2);
[HPBW2,SLL2,Nulls2,~] = Get_HPBW_SLL_Nulls(0,F_dB2);
%% �㷨�Ż�
SearchAgents_no=50; %��Ⱥ����
Max_iteration=200; %�趨����������
% ͬʱ���ڷ��Ⱥ���λ,x=[��λ������]
dim = 2*N;%ά��Ϊ4����x1-x4
ll = ones(1,N)*0.5;
ll(1,offarray) = 0;
lb = [zeros(1,N),ll];%�����±߽�
uu = ones(1,N)*2;
uu(1,offarray) = 0;
ub =[ones(1,N)*100,uu];%�����ϱ߽� 
fobj = @(x) funP(x,f,dx,dy,HPBW,SLL,Nulls);
%%  1.��ȸ
t1 = cputime;
[Best_score1,~,~]=SSA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %��ʼ�Ż�
t2 = cputime;
time_SSA = t2 - t1;

%% 2.Pso�㷨
t1 = cputime;
re_fitness_function=@(indivi)fitness_function(indivi,f,dx,dy,HPBW,SLL,Nulls);
PSOparams= [5 Max_iteration 50 3 2 1 0.4 1500 1e-25 1000 NaN 0 0];
range = [lb',ub']; %�����仯��Χ(��ɾ���)
Max_V = 0.5*(range(:,2)-range(:,1));  %����ٶ�ȡ�仯��Χ��10%~20%
[Best_score2]=pso_Trelea_vectorized(re_fitness_function,dim,Max_V,range,0,PSOparams);
Best_score2 = Best_score2';
t2 = cputime;
time_PSO = t2 - t1;
% figure(1)
% plot(SSA_curve,'Color','r','linewidth',1.5)
%% 3.�Ż���ȸ
t1 = cputime;
[Best_score3,~,~]=SSA_lv(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %��ʼ�Ż�
t2 = cputime;
time_MSSA = t2 - t1;

%% �ȶ��Ż�ǰ��Ч��
% SSA
Phase_SSA = Best_score1(1,1:N)';
Amp_SSA = Best_score1(1,N+1:end-1)';
[F_dB_SSA] = PlanerArrayPattern(f,Amp_SSA,Phase_SSA,dx,dy);
[HPBW_SSA,SLL_SSA,Nulls_SSA,~] = Get_HPBW_SLL_Nulls(0,F_dB_SSA);

% PSO
Phase_PSO = Best_score2(1,1:N)';
Amp_PSO = Best_score2(1,N+1:end-1)';
[F_dB_PSO] = PlanerArrayPattern(f,Amp_PSO,Phase_PSO,dx,dy);
[HPBW_PSO, SLL_PSO, Nulls_PSO,~] = Get_HPBW_SLL_Nulls(0,F_dB_PSO);

% MSSA
Phase_MSSA = Best_score3(1,1:N)';
Amp_MSSA = Best_score3(1,N+1:end-1)';
[F_dB_MSSA] = PlanerArrayPattern(f,Amp_MSSA,Phase_MSSA,dx,dy);
[HPBW_MSSA,SLL_MSSA,Nulls_MSSA,~] = Get_HPBW_SLL_Nulls(0,F_dB_MSSA);


figure(2)    %Phi = 0��ѿ�������
res = pi/360;
Phi = 0:res:2*pi-res;
Theta = 0:res:pi/2;
NumPhi = length(Phi);
plot(Theta*180/pi,F_dB(1,:),'k');
hold on
plot(-Theta*180/pi,F_dB(NumPhi/2+1,:),'k');
grid on

plot(Theta*180/pi,F_dB_SSA(1,:),'blue');
hold on
plot(-Theta*180/pi,F_dB_SSA(NumPhi/2+1,:),'blue');
hold on

plot(Theta*180/pi,F_dB_PSO(1,:),'red');
hold on
plot(-Theta*180/pi,F_dB_PSO(NumPhi/2+1,:),'red');
hold on

plot(Theta*180/pi,F_dB_MSSA(1,:),'green');
hold on
plot(-Theta*180/pi,F_dB_MSSA(NumPhi/2+1,:),'green');
hold on

hold off
xlim([-90,90]);
% ylim([-40,0]);
% title('Phi = 0��ѿ�������');
xlabel('Theta(deg)');
ylabel('F(dB)');

figure(10)
title('Ŀ�꺯�����');
plot(Xresult2,'red');
hold on
plot(Xresult,'blue');
hold on
plot(Xresult3,'green');
hold on
legend('PSO','SSA','MSSA');
hold off;

fprintf('time_SSA = %f\n',time_SSA);
fprintf('time_PSO = %f\n',time_PSO);
fprintf('time_MSSA = %f\n',time_MSSA);
