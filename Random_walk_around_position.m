% 随机游走策略
function [RWs]=Random_walk_around_position(Dim,max_iter,lb, ub,position,current_iter)


I=1; % 系数
%系数根据迭代次数进行改变，迭代次数越大，搜索范围越小，寻优更精细
if current_iter>max_iter/10
    I=1+100*(current_iter/max_iter);
end

if current_iter>max_iter/2
    I=1+1000*(current_iter/max_iter);
end

if current_iter>max_iter*(3/4)
    I=1+10000*(current_iter/max_iter);
end

if current_iter>max_iter*(0.9)
    I=1+100000*(current_iter/max_iter);
end

if current_iter>max_iter*(0.95)
    I=1+1000000*(current_iter/max_iter);
end


% 根据系数改变搜索范围，迭代次数越大，范围越小。
lb=lb/(I);
ub=ub/(I); 

% 将范围移动到目标附近
if rand<0.5
    lb=lb+position; % Equation (2.8) in the paper
else
    lb=-lb+position;
end

if rand>=0.5
    ub=ub+position; 
else
    ub=-ub+position;
end

% 随机游走
for i=1:Dim
    X = [0 cumsum(2*(rand(max_iter,1)>0.5)-1)']; 
    %[a b]--->[c d]
    a=min(X);
    b=max(X);
    c=lb(i);
    d=ub(i);      
    X_norm=((X-a).*(d-c))./(b-a)+c; 
    RWs(:,i)=X_norm;
end

