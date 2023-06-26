%2023年5月26日
%作者：3200105710
%DMC控制算法实现与仿真

%初始化，清空变量及程序区
clc
clear all

% 建立被控对象传递函数，如传递函数为G(s)=1/(50s^2+15s+1)
num=[8611.77];%分子
den=[1 1.6 274.1 279.1 8611.7];%分母
G=tf(num,den);%创建传递函数模型

%定义并初始化超参数及状态变量，通过仿真调优确定
P=20; % 优化时域
M=1; % 控制时域
N=40; % 截断时域,此处需满足PPT49页中M<=P<=N的范围要求
Ts=2; % 根据对象特性，选择合适的采样周期
steps=100; % 仿真长度

% 根据传递函数得到离散状态空间方程
[As,Bs,Cs,Ds]=tf2ss(num,den);
[Ad,Bd]=c2d(As,Bs,Ts);
Xs0=[0 0 0 0]';

% 利用step函数获取传递函数阶跃响应，即求取PPT50页中DMC动态矩阵内部各参数
[a0,t]=step(G,0:Ts:(N-1)*Ts);

% 构造动态矩阵A，矩阵A维数即为P*M，优化时域*控制时域
As=zeros(P,M);%初始化
%A=[a1;a2 a1;a3 a2 a1;......;am am-1 ... a1;.......;ap ap-1 ...ap-m+1]p*m
As(:,1)=a0(1:P);%填充第一列
%观察每一列均是前一列向下平移单位1，向上填充0的结果，用循环实现该过程
for i=1:P
    for j=2:M
        if i>=j
            As(i,j)=As(i-1,j-1);
        end
    end
end

% 构造误差校正向量
h=0.5*ones(N,1);
h(1)=1;
s=zeros(N,N);

%构造转移矩阵
for i=1:N-1
    s(i,i+1)=1;
end
s(N,N)=1;


% 离线计算向量d，公式为PPT58页d^T=c^T(A^TQA+R)^-1A^TQ
Q=eye(P);
R=0*eye(M);
c=zeros(M,1);
c(1)=1;
d=c'*(As'*Q*As+R)^-1*As'*Q;

%在线计算部分变量初始化
yr=ones(P,1); % 参考轨迹
y0=zeros(N,1); % 模型预测值
y=zeros(steps,1); % 实际输出值
u=zeros(steps,1); % 系统控制量

% 在线计算首步计算，具体原理与见下方部分计算
xs1=Ad*Xs0;
y(1)=Cs*xs1;
Xs0=xs1;
ycor=y0+h*(y(1)-y0(1));
y0=s*ycor;
du=d*(yr-y0(1:P));
y0=y0+a0*du;
u(1)=du;

% 在线计算滚动优化
for k=2:steps
    %根据状态空间模型计算T时刻xk和yk作为测量到的实际输出值
    xs1=Ad*Xs0+Bd*u(k-1);
    y(k)=Cs*xs1+Ds*u(k-1);
    Xs0=xs1;
    ycor=y0+h*(y(k)-y0(1));%误差修正后的输出预测值，公式为YN1cor(k+1)=YN1(k)+h*e(k+1)
    y0=s*ycor;%计算出预测值，公式YN0(k+1)=S*YNcor(k+1)
    du=d*(yr-y0(1:P));%对应公式du(k)=c^TU(k)=d^T(W(k)-Y(k))
    y0=y0+a0*du;%模型预测公式2，计算出YN1(k)=YN0(k)+adu(k)
    u(k)=u(k-1)+du;%更新u(k)
end

% 绘制图形
figure(1);
subplot(2,1,1);
plot(y,'linewidth',2);
title('Outputs');
xlabel('t');
ylabel('y');
hold on
grid on;
subplot(2,1,2)
plot(u,'linewidth',2);
title('Manipulator Variables');
xlabel('t');
ylabel('u');
grid on;