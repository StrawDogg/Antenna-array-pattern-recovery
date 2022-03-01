figure(10)
Xresult =[1, 1, 1, 1, 1];
Xresult1 =[1, 1, 1, 1, 1];
Xresult2 =[1, 1, 1, 1, 1];
Xresult3 =[1, 1, 1, 1, 1];
x = [1 2 3 4 5];
plot(x,Xresult,'k');
hold on
plot(x,Xresult1,'red');
hold on
plot(x,Xresult2,'blue');
hold on
plot(x,Xresult3,'green');
hold on
legend('Uniform','PSO','SSA','MSSA');