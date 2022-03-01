%结合PlanerArrayPattern获取指定Phi角度Phi0对应的平面内的HPBWW,SLL和所有零点
%Phi0取值范围为 0:res:180-res ,因为Phi0和Phi0+180是同一个平面所以上限不是360-res
%res_Ipt为插值后的角度分辨率,因为若res太小，方向图计算过于耗时，采用大res计算方向图再插值的办法

%HPBW为一个数
%SLL为一个数
%Nulls为一个两行的矩阵，第一行是所有零点的角度 ，第二行是所有零点的深度，即功率,当前的插值方法对提高零点精度似乎无效

function [HPBW,SLL,Nulls,F_dB_Nmlz_Phi0_Ipt] = Get_HPBW_SLL_Nulls(Phi0,F_dB_Nmlz,res_Ipt)

if nargin == 2
    res_Ipt = 0.1;
end

%获取Phi0平面的方向图
[NumPhi,~] = size(F_dB_Nmlz);
res = 360/NumPhi;           %角度分辨率
N1 = Phi0/res + 1;          %Phi0对应的F_dB_Nmlz中的行数
N2 = (180+Phi0)/res + 1;    %Phi0+180对应的F_dB_Nmlz中的行数
Theta = -90:res:90;         %俯仰角，负的俯仰角对应Phi0+180的半平面
F_dB_Nmlz_Phi0 = [F_dB_Nmlz(N2,end:-1:2),F_dB_Nmlz(N1,:)];  %Phi0对应的平面上的波束功率，其元素与Theta向量的元素一一对应

%为提高分辨率，对F_dB_Nmlz_Phi0插值
Theta_Ipt = -90:res_Ipt:90;
% Theta_Ipt = setdiff(Theta_Ipt,Theta);
try
 F_dB_Nmlz_Phi0_Ipt = interp1(Theta,F_dB_Nmlz_Phi0,Theta_Ipt,'pchip'); %对F_dB_Nmlz_Phi0插值
catch ME
    rethrow(ME);
end

%获取HPBW
N_Fmax = find(F_dB_Nmlz_Phi0_Ipt == max(F_dB_Nmlz_Phi0_Ipt));  %主瓣最高功率点对应的序号
%容易找到多个点，所以取个平均__
N_Fmax =N_Fmax(1,1);
N_HBPW_L = N_Fmax;           %主瓣左侧的半功率点对应的序号
% N_HBPW_L
N_HBPW_R = N_Fmax;           %主瓣右侧的半功率点对应的序号
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
%获取SLL
Peaks = findpeaks(F_dB_Nmlz_Phi0_Ipt); %找到主瓣和所有旁瓣的最高功率点
[Peaks,~] = sort(Peaks);          %从小到大排序
SLL = Peaks(end) - Peaks(end-1);  %用主瓣功率减去最大的旁瓣功率

%获取零点
Nulls_Depth = -findpeaks(-F_dB_Nmlz_Phi0_Ipt); %找到所有最低功率点，即零点

% 筛选零点
shaixuan = find(Nulls_Depth < -40);
Nulls_Depth = Nulls_Depth(shaixuan);

Nulls_index = zeros(1,length(Nulls_Depth));
Nulls_Angle = zeros(1,length(Nulls_Depth));   %初始化零点角度向量
for i = 1:length(Nulls_Depth)
    index = find(F_dB_Nmlz_Phi0_Ipt==Nulls_Depth(i));  %找到零点的角度
    Nulls_index(i) = index(1);
    Nulls_Angle(i) = Theta_Ipt(Nulls_index(i));
end
Nulls = [Nulls_index;Nulls_Depth;Nulls_Angle];


%绘图
% figure
% plot(Theta_Ipt,F_dB_Nmlz_Phi0_Ipt);
% xlim([-90,90]);
% ylim([-40,0]);
end
