%效果验证
close all
clearvars -except Theta   %只保留辨识出来的参数进行检验

%% 取参数的最后五个数据的平均值作为检验参数
sizeiden=size(Theta);
ls=Theta(1:4,sizeiden(2)-4:sizeiden(2));
ls=mean(ls,2);

%% 读取数据与处理
load rd_2.mat
datalength=data(1,1); %数据总长
delay=18;  %观察所得延时时间 test1时延为20，test2时延为8

length=datalength-delay-2;  %真实输出数据的长度
u=data(2:datalength+1,2);%输入
y=data(2:datalength+1,3); %输出
u0=u(3:length+2);%u(t)
u1=u(2:length+1);%u(t-1)
u2=u(1:length);%u(t-2)

%% 对原始数据进行拟合，得到辨识后的输出
y2(1:delay+2,1)=y(1:delay+2);
for k=3:datalength
    i=k-delay-2;  %保证从延时之后的时间开始拟合
    if(i>0)       
        y2(k,1)=[-y2(k-1,1) -y2(k-2,1) u1(i) u2(i)]*ls;
    else          %延时之前的输入设定为初始输入
        y2(k,1)=[-y2(k-1,1) -y2(k-2,1) 6.5 6.5]*ls;   
    end
end
%% 原辨识系统的结果对比
t=1:datalength;
figure(1);
plot(t,y);hold on;
plot(t,y2);hold off;
grid;
title("渐消记忆递推验证对比");
legend('实际输出','辨识系统的输出');

figure(2)
yerror=y-y2;    %计算辨识曲线与原曲线的误差
plot(t,yerror);grid; J=yerror'*yerror;%检验误差J=e^2
Jvc=J/datalength;   
title("渐消记忆递推辨识系统与实际系统误差曲线 e");    
fprintf("检验误差：J = %.10f\n",J);
fprintf("检验方差：Jvc = %.10f\n",Jvc);
