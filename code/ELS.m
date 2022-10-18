%增广最小二乘递推算法JL21010001 张笑一
clear;clc;close all
%% 导入数据
load train.mat
datalength=data(1,1);%数据总长
u=data(2:datalength,2);%输入数据
v=randn(1,datalength);%构造随机扰动
delay=20;%时延
length=datalength-delay;%减去时延后总长
y=data(2+delay:datalength,3);

%% 递推求解
P=10^6*eye(6); %选取P的初始值为100-10^6之间，一般选10^6
Theta=zeros(6,length-2);     %参数的估计值，存放中间过程估值
% K=[10;10;10;10;10;10];     %初始化K矩阵
for i=3:length-1
    h=[-y(i-1);-y(i-2);u(i-1);u(i-2);v(i-1);v(i-2)];%将二阶系统每一次的输入、输出、扰动导入进行迭代
    K=P*h*inv(h'*P*h+1);    %更新K值
    P=(eye(6)-K*h')*P;      %更新P值
    Theta(:,i-1)=Theta(:,i-2)+K*(y(i)-h'*Theta(: ,i-2));%计算参数向量[a1,a2,b1,b2,c1,c2]，全部放在Theta矩阵中
    v(i)=y(i)-h'*Theta(:,i-1);  %更新扰动
end

%% 计算误差和绘制参数辨识结果
L=length-2;%真值数据长度
y1=y(2:length-1);%真实值
y2=zeros(L,1);%估计值
a=size(y2);
for i=1:2 %1-2
 y2=y2-Theta(i,a(1))*y(3-i:length-i);%先计算输出的前两阶
end 
for i=1:2 
 y2=y2+Theta(i+2,a(1))*u(3-i:length-i);%计算输入的前两阶
end
error=y1-y2;%残差向量
J=error'*error;
J2=J/L;%残差方差
fprintf("增广最小二乘法\n");
fprintf("误差：J=%.10f\n",J);
fprintf("方差：J2=%.10f\n",J2);

%绘制辨识参数曲线
sizeiden=size(Theta);
i=1:sizeiden(2);
figure (1)

plot(i, Theta(1, :), i,Theta(2,:) , i,Theta(3,:), i, Theta(4,:), i, Theta(5,:), i,Theta(6,:))
title('增广最小二乘法')
legend('a1','a2','b1','b2','c1','c2');



