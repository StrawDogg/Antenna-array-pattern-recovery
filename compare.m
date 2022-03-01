clc;clear;
close all;
clear funP
import Config.*

load w
%% ���߻�������
f = 3.5e9; % Ƶ��
N = Config.N; % ��Ԫ����
dx = 3e8/f/2; % ��Ԫ���Ϊ������һ��
dy = 0.043;  %��Ϊy����ֻ��1�����ߵ�Ԫ����ֵ����ν

%% ��������λ��
Amp = ones(N,1);
Phase =  randn(N,1)*10 + 90; % ��̬
[F_dB] = PlanerArrayPattern(f,Amp,Phase,dx,dy);
[HPBW,SLL,Nulls,~] = Get_HPBW_SLL_Nulls(0,F_dB);
% %% ����ʧЧ
% offarray=[1 2 18 19]
% 
% Amp2 = ones(N,1);
% Amp2(offarray,1) = 0;
% [F_dB2] = PlanerArrayPattern(f,Amp2,Phase,dx,dy);
[HPBW2,SLL2,Nulls2,~] = Get_HPBW_SLL_Nulls(0,F_dB);

Amp2 = ones(N,1);
Amp2 = abs(w*Amp2);
[F_dB2] = PlanerArrayPattern(f,Amp2,Phase,dx,dy);

figure(2)    %Phi = 0��ѿ�������
res = pi/360;
Phi = 0:res:2*pi-res;
Theta = 0:res:pi/2;
NumPhi = length(Phi);
plot(Theta*180/pi,F_dB(1,:),'blue');
hold on
plot(-Theta*180/pi,F_dB(NumPhi/2+1,:),'blue');
hold on
plot(Theta*180/pi,F_dB2(1,:),'red');
hold on
plot(-Theta*180/pi,F_dB2(NumPhi/2+1,:),'red');
hold off

 