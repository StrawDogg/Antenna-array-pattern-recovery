%1*8������x����1��Ԫ��y����8��Ԫ�����Ϊ0.043m�����Ⱦ�Ϊ1V��������,��λ��Ϊ90deg�����������3.5GHz�ķ���ͼ
clc,clear;
f = 3.5e9;
Amp = ones(8,1);
Phase = [0:7]'*10;
Phase = zeros(8,1);
dx = 0.043;
dy = 0.043;  %��Ϊy����ֻ��1�����ߵ�Ԫ����ֵ����ν
[F_dB] = PlanerArrayPattern(f,Amp,Phase,dx,dy);
[HPBW,SLL,Nulls] = Get_HPBW_SLL_Nulls(0,F_dB);