% ������߲���
function [RWs]=Random_walk_around_position(Dim,max_iter,lb, ub,position,current_iter)


I=1; % ϵ��
%ϵ�����ݵ����������иı䣬��������Խ��������ΧԽС��Ѱ�Ÿ���ϸ
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


% ����ϵ���ı�������Χ����������Խ�󣬷�ΧԽС��
lb=lb/(I);
ub=ub/(I); 

% ����Χ�ƶ���Ŀ�긽��
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

% �������
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

