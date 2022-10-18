%渐消记忆递推最小二乘JL21060004 廖锦涛
clc,clear;close all;

load train.mat
datalength=data(1,1);%数据总长
delay=20;%时延
mu=0.996;%遗忘因子
length=datalength-delay;%减去时延后的数据总长
u=data(2:length+1,2);%输入
y=data(2+delay:datalength+1,3);%输出

%递推算法求解参数
%参数初始化
theta = zeros(4,1);
Theta = zeros(4,1);
J=0;
P = 10^6*eye(4);%P一般选取10^6
for k = 3:length-1
    h = [-y(k-1) -y(k-2) u(k-1) u(k-2)]';%将二阶系统每一次的输入、输出导入进行迭代
    K = P*h/(mu + h'*P*h);      %更新K
    P = (P - K*h'*P)/mu;        %更新P
    J=J+(y(k) - h'*theta)^2/(mu + h'*P*h);  %计算准则函数
    theta = theta + K*(y(k) - h'*theta);    %更新辨识参数a1,a2,b1,b2
    Theta = [Theta,theta];                  %存放辨识参数
end

%%计算误差与绘制参数曲线
L=length-2;%真值数据长度
y1=y(3:length);%真实值
y2=zeros(L,1);%估计值
a=size(y2);
for i=1:2 
 y2=y2-Theta(i,a(1))*y(3-i:length-i);%先计算输出的辨识树池
end
for i=1:2 
 y2=y2+Theta(i+2,a(1))*u(3-i:length-i);%计算输入的辨识输出
end
error=y1-y2;%残差向量
J=error'*error;
J2=J/L;%残差方差
if(mu==1)
     fprintf("递推最小二乘法 RLS\n");
else fprintf("渐消记忆递推最小二乘法 MLS\n");
end
fprintf("误差：J = %.10f\n",J);
fprintf("方差：J2 = %.10f\n",J2);

sizeiden=size(Theta);
k=1:sizeiden(2);
plot(k,Theta');grid;
if(mu==1)
     title("加扰动递推最小二乘参数变化 ");
else title('加扰动渐消记忆递推最小二乘法')
end
legend('a1','a2','b1','b2');







