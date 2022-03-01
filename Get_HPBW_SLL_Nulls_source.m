%���PlanerArrayPattern��ȡָ��Phi�Ƕ�Phi0��Ӧ��ƽ���ڵ�HPBWW,SLL���������
%Phi0ȡֵ��ΧΪ 0:res:180-res ,��ΪPhi0��Phi0+180��ͬһ��ƽ���������޲���360-res
%res_IptΪ��ֵ��ĽǶȷֱ���,��Ϊ��res̫С������ͼ������ں�ʱ�����ô�res���㷽��ͼ�ٲ�ֵ�İ취

%HPBWΪһ����
%SLLΪһ����
%NullsΪһ�����еľ��󣬵�һ�����������ĽǶ� ���ڶ���������������ȣ�������,��ǰ�Ĳ�ֵ�����������㾫���ƺ���Ч

function [HPBW,SLL,Nulls,F_dB_Nmlz_Phi0_Ipt] = Get_HPBW_SLL_Nulls(Phi0,F_dB_Nmlz,res_Ipt)

if nargin == 2
    res_Ipt = 0.1;
end

%��ȡPhi0ƽ��ķ���ͼ
[NumPhi,~] = size(F_dB_Nmlz);
res = 360/NumPhi;           %�Ƕȷֱ���
N1 = Phi0/res + 1;          %Phi0��Ӧ��F_dB_Nmlz�е�����
N2 = (180+Phi0)/res + 1;    %Phi0+180��Ӧ��F_dB_Nmlz�е�����
Theta = -90:res:90;         %�����ǣ����ĸ����Ƕ�ӦPhi0+180�İ�ƽ��
F_dB_Nmlz_Phi0 = [F_dB_Nmlz(N2,end:-1:2),F_dB_Nmlz(N1,:)];  %Phi0��Ӧ��ƽ���ϵĲ������ʣ���Ԫ����Theta������Ԫ��һһ��Ӧ

%Ϊ��߷ֱ��ʣ���F_dB_Nmlz_Phi0��ֵ
Theta_Ipt = -90:res_Ipt:90;
% Theta_Ipt = setdiff(Theta_Ipt,Theta);
try
 F_dB_Nmlz_Phi0_Ipt = interp1(Theta,F_dB_Nmlz_Phi0,Theta_Ipt,'pchip'); %��F_dB_Nmlz_Phi0��ֵ
catch ME
    rethrow(ME);
end

%��ȡHPBW
N_Fmax = find(F_dB_Nmlz_Phi0_Ipt == max(F_dB_Nmlz_Phi0_Ipt));  %������߹��ʵ��Ӧ�����
%�����ҵ�����㣬����ȡ��ƽ��__
N_Fmax =N_Fmax(1,1);
N_HBPW_L = N_Fmax;           %�������İ빦�ʵ��Ӧ�����
% N_HBPW_L
N_HBPW_R = N_Fmax;           %�����Ҳ�İ빦�ʵ��Ӧ�����
[~,iptSize] = size(F_dB_Nmlz_Phi0_Ipt);
while (F_dB_Nmlz_Phi0_Ipt(N_HBPW_L) > -3)
    N_HBPW_L = N_HBPW_L - 1;
    if N_HBPW_L < 1
        N_HBPW_L = iptSize;
    end
end
while (F_dB_Nmlz_Phi0_Ipt(N_HBPW_R) > -3)
    N_HBPW_R = N_HBPW_R + 1;
    if N_HBPW_R > iptSize
        N_HBPW_R = 1;
    end
end    
HPBW = Theta_Ipt(N_HBPW_R) - Theta_Ipt(N_HBPW_L);
if HPBW < 0
    HPBW = HPBW + 180;
end
%��ȡSLL
Peaks = findpeaks(F_dB_Nmlz_Phi0_Ipt); %�ҵ�����������԰����߹��ʵ�
[Peaks,~] = sort(Peaks);          %��С��������
SLL = Peaks(end) - Peaks(end-1);  %�����깦�ʼ�ȥ�����԰깦��

%��ȡ���
Nulls_Depth = -findpeaks(-F_dB_Nmlz_Phi0_Ipt); %�ҵ�������͹��ʵ㣬�����

% ɸѡ���
shaixuan = find(Nulls_Depth < -40);
Nulls_Depth = Nulls_Depth(shaixuan);

Nulls_index = zeros(1,length(Nulls_Depth));
Nulls_Angle = zeros(1,length(Nulls_Depth));   %��ʼ�����Ƕ�����
for i = 1:length(Nulls_Depth)
    index = find(F_dB_Nmlz_Phi0_Ipt==Nulls_Depth(i));  %�ҵ����ĽǶ�
    Nulls_index(i) = index(1);
    Nulls_Angle(i) = Theta_Ipt(Nulls_index(i));
end
Nulls = [Nulls_index;Nulls_Depth;Nulls_Angle];


%��ͼ
% figure
% plot(Theta_Ipt,F_dB_Nmlz_Phi0_Ipt);
% xlim([-90,90]);
% ylim([-40,0]);
end
