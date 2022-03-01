%1*8天线阵，x方向1单元，y方向8单元，间距为0.043m，幅度均为1V带随机误差,相位差为90deg带随机误差，绘制3.5GHz的方向图
clc,clear;
f = 3.5e9;
Amp = ones(8,1);
Phase = [0:7]'*10;
Phase = zeros(8,1);
dx = 0.043;
dy = 0.043;  %因为y方向只有1个天线单元，此值无所谓
[F_dB] = PlanerArrayPattern(f,Amp,Phase,dx,dy);
[HPBW,SLL,Nulls] = Get_HPBW_SLL_Nulls(0,F_dB);