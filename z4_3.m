clc;clear;
close all;
clear funP
import Config.*

global Xresult;
Xresult = [];
%% 天线基础设置
f = 3.5e9; % 频率
N = Config.N; % 单元个数
dx = 3e8/f/2; % 单元间距为波长的一半
dy = 0.043;  %因为y方向只有1个天线单元，此值无所谓
%% 激励和相位差
Amp = ones(N,1);
load w
Amp = abs(w*Amp);
Phase = [0:N-1]'*10;
Phase = ones(N,1) *90;
Phase =  randn(N,1)*10 + 90; % 正态
[F_dB] = PlanerArrayPattern(f,Amp,Phase,dx,dy);
% [HPBW,SLL] = Get_HPBW_and_SLL(0,F_dB);
[HPBW,SLL,Nulls,~] = Get_HPBW_SLL_Nulls(0,F_dB);
%% 设置单元失效
Amp2 = Amp;
offarray=[10 18];
Amp2(offarray,1) = 0;
[F_dB2] = PlanerArrayPattern(f,Amp2,Phase,dx,dy);
% [HPBW2,SLL2] = Get_HPBW_and_SLL(0,F_dB2);
[HPBW2,SLL2,Nulls2,~] = Get_HPBW_SLL_Nulls(0,F_dB);
%% 麻雀算法优化
SearchAgents_no=100; %种群数量
Max_iteration=300; %设定最大迭代次数
% 同时调节幅度和相位,x=[相位，幅度]
dim = 2*N;%维度为4，即x1-x4
ll = ones(1,N)*0.01;
ll(1,offarray) = 0;
lb = [zeros(1,N),ll];%参数下边界
uu = ones(1,N)*4;
uu(1,offarray) = 0;
ub =[ones(1,N)*360,uu];%参数上边界
fobj = @(x) funP(x,f,dx,dy,HPBW,SLL,Nulls);
[Best_score,Best_pos,SSA_curve]=SSA_lv(SearchAgents_no,Max_iteration,lb,ub,dim,fobj); %开始优化
figure(1)
plot(SSA_curve,'Color','r','linewidth',1.5)
%% 比对优化前后效果
Phase = Best_score(1,1:N)';
Amp2 = Best_score(1,N+1:end)';
[F_dB3] = PlanerArrayPattern(f,Amp2,Phase,dx,dy);
[HPBW3,SLL3,Nulls3,~] = Get_HPBW_SLL_Nulls(0,F_dB);
figure(2)    %Phi = 0°笛卡尔坐标
res = pi/360;
Phi = 0:res:2*pi-res;
Theta = 0:res:pi/2;
NumPhi = length(Phi);

plot(Theta*180/pi,F_dB(1,:),'blue');
hold on
plot(-Theta*180/pi,F_dB(NumPhi/2+1,:),'blue');
grid on

plot(Theta*180/pi,F_dB2(1,:),'red');
hold on
plot(-Theta*180/pi,F_dB2(NumPhi/2+1,:),'red');
hold on

plot(Theta*180/pi,F_dB3(1,:),'g');
hold on
plot(-Theta*180/pi,F_dB3(NumPhi/2+1,:),'g');
hold off

xlim([-90,90]);
% ylim([-40,0]);
% title('Phi = 0°笛卡尔坐标');
xlabel('Theta(deg)');
ylabel('F(dB)');

% 
% figure(3)
% title('目标函数结果');
% plot(Xresult);
% hold on