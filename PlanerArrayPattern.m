%根据xy平面天线阵的各单元的幅相，以及各天线单元间距，画出该天线阵的方向图（z>0)，并给出最大辐射方向

%DPhi:最大辐射方向的Phi角度,即方位角，单位为deg，取值范围[-180，180],+x为Phi=0deg,+y为Phi=90deg,
%DTheta:最大辐射方向的Theta角度,即俯仰角，单位为deg，取值范围[0，180]，+z为Theta=0deg，xy平面为Theta=90deg
%f:频率,单位Hz
%Amp:二维矩阵，每个单元的馈电幅度（V），若为行向量则表示y方向的一维线阵，若为列向量则表示x方向的一维线阵
%Phase:二维矩阵，每个单元的馈电相位，单位为deg,行列数和Amp相同
%dx:x方向每个单元间隔距离，单位为米
%dy:y方向每个单元间隔距离，单位为米
%res:求解的角度分辨率，缺省为0.5°

function [F_dB_Nmlz] = PlanerArrayPattern(f,Amp,Phase,dx,dy,res)

if nargin == 5
    res = pi/360;
end

[Numx,Numy] = size(Amp);
Phase = Phase/180*pi;
k = 2*pi*f/3e8;
Phi = 0:res:2*pi-res;
Theta = 0:res:pi/2;
NumPhi = length(Phi);
NumTheta = length(Theta);
F = zeros(NumPhi,NumTheta); %每个方向的辐射强度
rA = zeros(Numx,Numy);
for nx = 1:Numx
    for ny = 1:Numy
        rA(nx,ny) = sqrt(((nx-1)*dx)^2 + ((ny-1)*dy)^2);
    end
end
PhiA = zeros(Numx,Numy);
for nx = 1:Numx
    for ny = 1:Numy
        if nx==1 && ny==1
            PhiA(nx,ny) = 0;
        else
            PhiA(nx,ny) = acos( (nx-1)*dx / rA(nx,ny));
        end
    end
end
for i = [1,NumPhi/2+1]
% for i = 1:NumPhi
    for j = 1:NumTheta
        for nx = 1:Numx
            for ny = 1:Numy
                F(i,j) = F(i,j) + Amp(nx,ny)*exp(1j*(k*rA(nx,ny)*sin(Theta(j))*cos(Phi(i) - PhiA(nx,ny)) + Phase(nx,ny)));
            end
        end
    end
end
F = abs(F)/(Numx*Numy);    %归一化
F = F.^2;    %换算为功率
F_dB = 10*log10(F); %换算为dB
F_dB_Nmlz = F_dB - max(max(F_dB));  %归一化

%%
%求DPhi和DTheta
% if size(Amp,1) == 1        %y方向的线阵。线阵的最大辐射方向无穷多，取线阵所在的与xoy垂直的平面上的最大辐射方向
%     F90 = ([F(NumPhi*2/8+1,:),F(NumPhi*6/8+1,:)]);  %90°面上的辐射强度
%     F90Max = max(F90);
%     IndexF90Max = find(F90 == F90Max);
%     F90MaxNum = length(IndexF90Max);
%     DPhi = zeros(1,F90MaxNum);
%     DTheta = zeros(1,F90MaxNum);
%     for i = 1:F90MaxNum
%         if IndexF90Max(i) == 1 || IndexF90Max(i) == 1+NumTheta %辐射方向为+z
%             DPhi(i) = 0;
%             DTheta(i) = 0;
%         else
%             if IndexF90Max(i) >=2 && IndexF90Max(i) <= NumTheta
%                 DPhi(i) = 90;
%                 DTheta(i) = Theta(IndexF90Max)/pi*180;
%             else
%                 DPhi(i) = 270;
%                 DTheta(i) = Theta(IndexF90Max - NumTheta)/pi*180;
%             end
%         end
%     end
% else
%     if size(Amp,2) == 1                %x方向的线阵
%         F180 = ([F(NumPhi*0/8+1,:),F(NumPhi*4/8+1,:)]);  %0°面上的辐射强度
%         F180Max = max(F180);
%         IndexF180Max = find(F180 == F180Max);
%         F180MaxNum = length(IndexF180Max);
%         DPhi = zeros(1,F180MaxNum);
%         DTheta = zeros(1,F180MaxNum);
%         for i = 1:F180MaxNum
%             if IndexF180Max(i) == 1 || IndexF180Max(i) == 1+NumTheta %辐射方向为+z
%                 DPhi(i) = 0;
%                 DTheta(i) = 0;
%             else
%                 if IndexF180Max(i) >=2 && IndexF180Max(i) <= NumTheta
%                     DPhi = 0;
%                     DTheta = Theta(IndexF180Max(i))/pi*180;
%                 else
%                     DPhi = 180;
%                     DTheta = Theta(IndexF180Max(i) - NumTheta)/pi*180;
%                 end
%             end
%         end
%     else      %面阵
%         [IndexMaxF1,IndexMaxF2] = find(F == max(max(F)));
%         DPhi = Phi(IndexMaxF1)/pi*180;
%         DTheta = Theta(IndexMaxF2)/pi*180;
%     end        
% end

%%
%绘图

% figure(1)    %Phi = 0°极坐标
% polarplot(Theta,F_dB_Nmlz(1,:),'blue');
% hold on
% polarplot(-Theta,F_dB_Nmlz(NumPhi/2+1,:),'blue');
% hold off
% title('Phi = 0°极坐标');
% ax = gca;
% ax.ThetaZeroLocation = 'top';
% ax.ThetaDir = 'clockwise';
% ax.RLim = [-40,0];
% ax.ThetaLim = [-90,90];
% 
figure(1)    %Phi = 0°笛卡尔坐标
plot(Theta*180/pi,F_dB_Nmlz(1,:),'blue');
hold on
plot(-Theta*180/pi,F_dB_Nmlz(NumPhi/2+1,:),'blue');
grid on
hold off
xlim([-90,90]);
% ylim([-40,0]);
title('Phi = 0°笛卡尔坐标');
xlabel('Theta(deg)');
ylabel('归一化幅度(dB)');
% 
% figure(3)    %Phi = 90°极坐标
% polarplot(Theta,F_dB_Nmlz(NumPhi*2/8+1,:),'blue');
% hold on
% polarplot(-Theta,F_dB_Nmlz(NumPhi*6/8+1,:),'blue');
% hold off
% title('Phi = 90°极坐标');
% ax = gca;
% ax.ThetaZeroLocation = 'top';
% ax.ThetaDir = 'clockwise';
% ax.RLim = [-40,0];
% ax.ThetaLim = [-90,90];
% 
% figure(1)    %Phi = 90°笛卡尔坐标
% plot(Theta*180/pi,F_dB_Nmlz(NumPhi*2/8+1,:),'blue');
% hold on
% plot(-Theta*180/pi,F_dB_Nmlz(NumPhi*6/8+1,:),'blue');
% grid on
% hold off
% xlim([-90,90]);
% ylim([-40,0]);
% title('Phi = 0°笛卡尔坐标');
% xlabel('Theta(deg)');
% ylabel('归一化幅度(dB)');
% 
% Fx = zeros(NumPhi,NumTheta);
% Fy = zeros(NumPhi,NumTheta);
% Fz = zeros(NumPhi,NumTheta);
% for i = 1:NumPhi
%     for j = 1:NumTheta
%         [Fx(i,j),Fy(i,j),Fz(i,j)] = sph2cart(Phi(i),pi/2-Theta(j),F(i,j));
%     end
% end
% 
% figure(5)       %3D图
% mesh(Fx,Fy,Fz);
% xlabel('x');ylabel('y');zlabel('z');
% axis([-1,1,-1,1,0,1]);
% view([1,-1,1]);    %3D图的观察角度
% title('3D方向图');
% zlabel('归一化幅度（非dB）');