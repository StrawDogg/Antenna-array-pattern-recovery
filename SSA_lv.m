%_________________________________________________________________________%
% 麻雀优化算法             %
%_________________________________________________________________________%
function [Best_pos,Best_score,curve]=SSA_lv(pop,Max_iter,lb,ub,dim,fobj)

ST = 0.6;%预警值
PD = 0.5;%发现者的比列，剩下的是加入者
SD = 0.2;%意识到有危险麻雀的比重
global Xresult3;

PDNumber = round(pop*PD); %发现者数量
SDNumber = round(pop*SD);%意识到有危险麻雀数量
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%种群初始化
X0=initialization(pop,dim,ub,lb);
X = X0;
%计算初始适应度值
fitness = zeros(1,pop);
Cnt = 1;
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
   disp(Cnt);
   Cnt = Cnt + 1;
end
 [fitness, index]= sort(fitness);%排序
BestF = fitness(1);
WorstF = fitness(end);
GBestF = fitness(1);%全局最优适应度值
for i = 1:pop
    X(i,:) = X0(index(i),:);
end
curve=zeros(1,Max_iter);
GBestX = X(1,:);%全局最优位置
X_new = X;
for count = 1: Max_iter
    BestF = fitness(1);
    WorstF = fitness(end);
    R2 = rand(1);
   for j = 1:PDNumber
      if(R2<ST)
          X_new(j,:) = X(j,:).*exp(-j/(rand(1)*Max_iter));
      else
          X_new(j,:) = X(j,:) + randn()*ones(1,dim);
      end     
   end
   for j = PDNumber+1:pop
%        if(j>(pop/2))
        if(j>(pop - PDNumber)/2 + PDNumber)
          X_new(j,:)= randn().*exp((X(end,:) - X(j,:))/j^2);
       else
          %产生-1，1的随机数
          A = ones(1,dim);
          for a = 1:dim
            if(rand()>0.5)
                A(a) = -1;
            end
          end 
          AA = A'*inv(A*A');
          X_new(j,:)= X(1,:) + abs(X(j,:) - X(1,:)).*AA';
       end
   end
   Temp = randperm(pop);
   SDchooseIndex = Temp(1:SDNumber); 
   for j = 1:SDNumber
       if(fitness(SDchooseIndex(j))>BestF)
            % X_new(SDchooseIndex(j),:) = X(1,:) + randn().*abs(X(SDchooseIndex(j),:) - X(1,:));
           % 莱维
           lv = levyFlight();
           X_new(SDchooseIndex(j),:) = lv * X(1,:) + randn().*abs(X(SDchooseIndex(j),:) - lv * X(1,:));
       elseif(fitness(SDchooseIndex(j))== BestF)
           K = 2*rand() -1;
           X_new(SDchooseIndex(j),:) = X(SDchooseIndex(j),:) + K.*(abs( X(SDchooseIndex(j),:) - X(end,:))./(fitness(SDchooseIndex(j)) - fitness(end) + 10^-8));
       end
   end
 %边界控制
   for j = 1:pop 
       for a = 1: dim
           gamma = rand;
           if(X_new(j,a)>ub(a))
               X_new(j,a) = ub(a) - gamma * min(X_new(j,a) - ub(a), ub(a) - lb(a));
           end
           if(X_new(j,a)<lb(a))
               X_new(j,a) =lb(a) + gamma * min(lb(a) - X_new(j,a), ub(a) - lb(a));
           end
           if(X_new(j,a)>ub(a))
               X_new(j,a) = ub(a);
           end
       end
   end 
   %更新位置
   parfor j=1:pop
    fitness_new(j) = fobj(X_new(j,:));
   end
   for j = 1:pop
    if(fitness_new(j) < GBestF)
       GBestF = fitness_new(j);
        GBestX = X_new(j,:);   
    end
   end
   X = X_new;
   fitness = fitness_new;
    %排序更新
   [fitness, index]= sort(fitness);%排序
   BestF = fitness(1);
   Xresult3 = [Xresult3 BestF];
   WorstF = fitness(end);
   for j = 1:pop
      X(j,:) = X(index(j),:);
   end
   curve(i) = GBestF;
   fprintf('迭代次数：%d，全局最优值：%f\n',count,GBestF);
%    if (GBestF < 0.001)
%        break;
%    end
end
Best_pos =GBestX;
Best_score = curve(end);
end



