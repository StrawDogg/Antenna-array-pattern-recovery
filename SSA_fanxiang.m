%_________________________________________________________________________%
% ��ȸ�Ż��㷨             %
%_________________________________________________________________________%
function [Best_pos,Best_score,curve]=SSA(pop,Max_iter,lb,ub,dim,fobj)

ST = 0.6;%Ԥ��ֵ
PD = 0.5;%�����ߵı��У�ʣ�µ��Ǽ�����
SD = 0.2;%��ʶ����Σ����ȸ�ı���
global Xresult;

PDNumber = round(pop*PD); %����������
SDNumber = round(pop*SD);%��ʶ����Σ����ȸ����
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%��Ⱥ��ʼ��
X0=initialization(pop,dim,ub,lb);
 X = X0;
%�����ʼ��Ӧ��ֵ
fitness = zeros(1,pop);
Cnt = 1;
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
   disp(Cnt);
   Cnt = Cnt + 1;
end
 [fitness, index]= sort(fitness);%����
BestF = fitness(1);
WorstF = fitness(end);
GBestF = fitness(1);%ȫ��������Ӧ��ֵ
for i = 1:pop
    X(i,:) = X0(index(i),:);
end 
curve=zeros(1,Max_iter);
GBestX = X(1,:);%ȫ������λ��
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
          %����-1��1�������
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
           X_new(SDchooseIndex(j),:) = X(1,:) + randn().*abs(X(SDchooseIndex(j),:) - X(1,:));
       elseif(fitness(SDchooseIndex(j))== BestF)
           K = 2*rand() -1;
           X_new(SDchooseIndex(j),:) = X(SDchooseIndex(j),:) + K.*(abs( X(SDchooseIndex(j),:) - X(end,:))./(fitness(SDchooseIndex(j)) - fitness(end) + 10^-8));
       end
   end
   % ����ѧϰ����
   X_star = rand * (lb + ub) - X_new;
   X_new = [X_star; X_new];
   %�߽����
   for j = 1:pop * 2
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
   %����λ��
   for j=1:pop*2
    fitness_new(j) = fobj(X_new(j,:));
   end
%    for j = 1:pop*2
%     if(fitness_new(j) < GBestF)
%        GBestF = fitness_new(j); % ���ŵ���ʺ϶Ⱥ͸õ�ķ���
%        GBestX = X_new(j,:);   
%     end
%    end
   % X = X_new;
   % fitness = fitness_new;
    %�������
   [fitness_new, index]= sort(fitness_new);%����
   if(fitness_new(1) < GBestF)
       GBestF = fitness_new(1); % ���ŵ���ʺ϶Ⱥ͸õ�ķ���
       GBestX = X_new(index(1),:);   
   end
   
   BestF = fitness_new(1);
   BestX = X_new(index(1));
   Xresult = [Xresult GBestF];  % ��¼ȫ�����ŵ�
   
   % WorstF = fitness_new(end);
   % ����X
   rand_num = randsample(pop * 2,pop,false);
   [rand_num, ~] = sort(rand_num);
   for j = 1:pop
      X(j,:) = X_new(index(rand_num(j)),:);
      fitness(j) = fitness_new(rand_num(j));
   end
%    X(1,:) = BestX;
%    fitness(1) = BestF;
   
   X_new = X;
   curve(i) = GBestF;
   fprintf('����������%d��ȫ������ֵ��%f\n',count,GBestF);
%    if (GBestF < 0.001)
%        break;
%    end
end
Best_pos =GBestX;
Best_score = curve(end);
end



