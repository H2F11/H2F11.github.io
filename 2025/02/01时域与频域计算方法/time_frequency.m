clc;
clear all;
close all;

%% 结构参数
m = 120;       % 质量 (t) 
k = 3e5;       % 刚度 (N/mm)
w = sqrt(k/m); % 结构自振频率 (rad/s)
kesi = 0.05;   % 阻尼比
c = 2 * m * w * kesi;  % 阻尼系数 (N·s/mm)

%% 读取地震加速度数据
time_step = 0.005;           % 采样间隔 (s)
ground_acc_data = load('ElCentro.txt');  % 读取加速度数据
% 将地震加速度（假设单位为 g）转换为惯性力 (N)
% 注意：9.8 表示重力加速度 (m/s^2)，乘 1e3 是为了与 k (N/mm) 单位匹配
F_t = - m * ground_acc_data * 9.8 * 1e3;  % 地面作用力
N = length(F_t);             % 采样点数
time = (0:N-1) * time_step;  % 时间向量

%% 频域方法计算相对位移
% 计算 FFT（快速傅里叶变换）
F_Ft = fft(F_t);

% 构造对应的角频率向量（本质是因为Nyquist频率）
halfN = floor(N/2); % 统一处理偶数和奇数情况
omega = 2*pi*[0:halfN, -halfN:-1]'/(N*time_step);

% 计算频响函数 H(omega)
H_omega = 1 ./ (-m*omega.^2 + 1i*c*omega + k);

% 计算位移的频谱：X(omega) = H(omega)*F(omega)
X_fft = H_omega .* F_Ft;

% 逆 FFT 得到时域相对位移响应
x_freq = ifft(X_fft, 'symmetric');

% 绘制频域得到的相对位移频谱（仅显示正频率部分）
f = (0:halfN-1)'/(N*time_step);
figure;
plot(f, abs(X_fft(1:halfN)));
title('相对位移频谱');
xlabel('频率 (Hz)');
ylabel('幅值');
grid on;
saveas(gcf,"./结构相对位移频谱.png");

%% 使用 Newmark-β 法进行时域计算（平均加速度法：beta=1/4, gamma=1/2）
beta = 1/4;
gamma = 1/2;
dt = time_step;
T = N;  % 时间步数

% 初始化响应变量
x_newmark = zeros(T, 1);  % 位移 (mm)
v_newmark = zeros(T, 1);  % 速度 (mm/s)
a_newmark = zeros(T, 1);  % 加速度 (mm/s^2)

% 初始条件（假设初始位移和速度为零）
x_newmark(1) = 0;
v_newmark(1) = 0;
a_newmark(1) = (F_t(1) - c*v_newmark(1) - k*x_newmark(1)) / m;

% 预先计算有效刚度
K_eff = k + c*(gamma/(beta*dt)) + m/(beta*dt^2);

% 循环计算每个时间步的响应
for i = 1:T-1
    % Predictor（预测）步：
    x_pred = x_newmark(i) + dt*v_newmark(i) + dt^2/2*(1-2*beta)*a_newmark(i);
    v_pred = v_newmark(i) + dt*(1-gamma)*a_newmark(i);
    
    % 计算有效荷载（相当于“残差”）
    R_eff = F_t(i+1) - k*x_pred - c*v_pred;
    
    % 计算位移修正量 Δx
    delta_x = R_eff / K_eff;
    
    % 更新响应：位移、速度、加速度
    x_newmark(i+1) = x_pred + delta_x;
    v_newmark(i+1) = v_pred + gamma/(beta*dt)*delta_x;
    a_newmark(i+1) = delta_x/(beta*dt^2);
end

%% 绘图比较两种方法的结果
figure;
subplot(2,1,1);
plot(time, x_freq, 'b');
title('频域方法得到的相对位移时程曲线');
xlabel('时间 (s)');
ylabel('位移 (mm)');
grid on;

subplot(2,1,2);
plot(time, x_newmark, 'r');
title('Newmark-β 方法得到的相对位移时程曲线');
xlabel('时间 (s)');
ylabel('位移 (mm)');
grid on;
saveas(gcf, '时域与频域方法所得相对位移时程曲线对比.png');