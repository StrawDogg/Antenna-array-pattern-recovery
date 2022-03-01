%����xyƽ��������ĸ���Ԫ�ķ��࣬�Լ������ߵ�Ԫ��࣬������������ķ���ͼ��z>0)�������������䷽��

%DPhi:�����䷽���Phi�Ƕ�,����λ�ǣ���λΪdeg��ȡֵ��Χ[-180��180],+xΪPhi=0deg,+yΪPhi=90deg,
%DTheta:�����䷽���Theta�Ƕ�,�������ǣ���λΪdeg��ȡֵ��Χ[0��180]��+zΪTheta=0deg��xyƽ��ΪTheta=90deg
%f:Ƶ��,��λHz
%Amp:��ά����ÿ����Ԫ��������ȣ�V������Ϊ���������ʾy�����һά������Ϊ���������ʾx�����һά����
%Phase:��ά����ÿ����Ԫ��������λ����λΪdeg,��������Amp��ͬ
%dx:x����ÿ����Ԫ������룬��λΪ��
%dy:y����ÿ����Ԫ������룬��λΪ��
%res:���ĽǶȷֱ��ʣ�ȱʡΪ0.5��

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
F = zeros(NumPhi,NumTheta); %ÿ������ķ���ǿ��
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
F = abs(F)/(Numx*Numy);    %��һ��
F = F.^2;    %����Ϊ����
F_dB = 10*log10(F); %����ΪdB
F_dB_Nmlz = F_dB - max(max(F_dB));  %��һ��

%%
%��DPhi��DTheta
% if size(Amp,1) == 1        %y�������������������䷽������࣬ȡ�������ڵ���xoy��ֱ��ƽ���ϵ������䷽��
%     F90 = ([F(NumPhi*2/8+1,:),F(NumPhi*6/8+1,:)]);  %90�����ϵķ���ǿ��
%     F90Max = max(F90);
%     IndexF90Max = find(F90 == F90Max);
%     F90MaxNum = length(IndexF90Max);
%     DPhi = zeros(1,F90MaxNum);
%     DTheta = zeros(1,F90MaxNum);
%     for i = 1:F90MaxNum
%         if IndexF90Max(i) == 1 || IndexF90Max(i) == 1+NumTheta %���䷽��Ϊ+z
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
%     if size(Amp,2) == 1                %x���������
%         F180 = ([F(NumPhi*0/8+1,:),F(NumPhi*4/8+1,:)]);  %0�����ϵķ���ǿ��
%         F180Max = max(F180);
%         IndexF180Max = find(F180 == F180Max);
%         F180MaxNum = length(IndexF180Max);
%         DPhi = zeros(1,F180MaxNum);
%         DTheta = zeros(1,F180MaxNum);
%         for i = 1:F180MaxNum
%             if IndexF180Max(i) == 1 || IndexF180Max(i) == 1+NumTheta %���䷽��Ϊ+z
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
%     else      %����
%         [IndexMaxF1,IndexMaxF2] = find(F == max(max(F)));
%         DPhi = Phi(IndexMaxF1)/pi*180;
%         DTheta = Theta(IndexMaxF2)/pi*180;
%     end        
% end

%%
%��ͼ

% figure(1)    %Phi = 0�㼫����
% polarplot(Theta,F_dB_Nmlz(1,:),'blue');
% hold on
% polarplot(-Theta,F_dB_Nmlz(NumPhi/2+1,:),'blue');
% hold off
% title('Phi = 0�㼫����');
% ax = gca;
% ax.ThetaZeroLocation = 'top';
% ax.ThetaDir = 'clockwise';
% ax.RLim = [-40,0];
% ax.ThetaLim = [-90,90];
% 
figure(1)    %Phi = 0��ѿ�������
plot(Theta*180/pi,F_dB_Nmlz(1,:),'blue');
hold on
plot(-Theta*180/pi,F_dB_Nmlz(NumPhi/2+1,:),'blue');
grid on
hold off
xlim([-90,90]);
% ylim([-40,0]);
title('Phi = 0��ѿ�������');
xlabel('Theta(deg)');
ylabel('��һ������(dB)');
% 
% figure(3)    %Phi = 90�㼫����
% polarplot(Theta,F_dB_Nmlz(NumPhi*2/8+1,:),'blue');
% hold on
% polarplot(-Theta,F_dB_Nmlz(NumPhi*6/8+1,:),'blue');
% hold off
% title('Phi = 90�㼫����');
% ax = gca;
% ax.ThetaZeroLocation = 'top';
% ax.ThetaDir = 'clockwise';
% ax.RLim = [-40,0];
% ax.ThetaLim = [-90,90];
% 
% figure(1)    %Phi = 90��ѿ�������
% plot(Theta*180/pi,F_dB_Nmlz(NumPhi*2/8+1,:),'blue');
% hold on
% plot(-Theta*180/pi,F_dB_Nmlz(NumPhi*6/8+1,:),'blue');
% grid on
% hold off
% xlim([-90,90]);
% ylim([-40,0]);
% title('Phi = 0��ѿ�������');
% xlabel('Theta(deg)');
% ylabel('��һ������(dB)');
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
% figure(5)       %3Dͼ
% mesh(Fx,Fy,Fz);
% xlabel('x');ylabel('y');zlabel('z');
% axis([-1,1,-1,1,0,1]);
% view([1,-1,1]);    %3Dͼ�Ĺ۲�Ƕ�
% title('3D����ͼ');
% zlabel('��һ�����ȣ���dB��');