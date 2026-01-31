tic;
% % 启动并行池
% if isempty(gcp('nocreate'))
%     parpool('local'); % 启动默认并行池
% end
% 
% %% 参数设置
% % 系统参数
% h = 3;                     % 社区数量
% user_per_community = 50;   % 每个社区的用户数量
% total_users = h * user_per_community;
% F = 200;                    % 视频资源数量
% Per_m = 0.3;               % 每个社区选作缓存节点的用户比例
% iu_count_per_community = round(user_per_community * Per_m);
% nf = 100;                    % 每个视频的片段数
% chi = [1.5 , 2 , 5 , 10 , 12]; % 分辨率级别，单位:Mbps（bit）
% chi_star = chi(end);       % 最高分辨率编码比特率
% Delta_t = 1;               % 视频片段时长，单位:秒
% v_chi_star = chi_star * Delta_t / 8; % 最高分辨率视频片段大小，单位:MB（Byte)
% 
% % 区域和距离参数
% region_size = 1500;        % 正方形区域边长，单位:米
% sbs_coverage = 150;        % SBS覆盖半径，单位:米
% iu_coverage = 50;          % IU D2D覆盖半径，单位:米
% community_radius = 150;    % 社区半径，单位:米
% max_movement_dist = 2;    % 相邻时间点之间的最大移动距离
% 
% % 时间参数
% T_small = 100;              % 小时间尺度数量
% T_large = 1;               % 大时间尺度数量
% 
% % Zipf参数
% gamma_m = 0.56;             % Zipf分布参数：越接近1，最热门的视频被请求概率越大
% 
% % 社交感知图参数
% alpha = 0.2;               % 亲密度权重
% beta = 0.4;                % 偏好相似度权重
% gamma = 0.4;               % 预测重要程度权重
% 
% % 通信参数
% B = 5;                    % 带宽，单位:MHz
% P_sbs = 1;               % SBS发射功率，单位:W 30dbm
% P_iu = 0.3;               % IU发射功率，单位:W 23dbm
% N0 = 1e-7;                 % 噪声功率，单位:W
% K = 1;                     % 系统参数
% epsilon = 3;               % 路径损失指数
% slowfading_dB = 3;% 慢衰落（阴影衰落）的标准差 (dB)
% 
% %% QoE模型参数设置
% alpha_qoe = 0.45;    % 视频分辨率权重
% beta_qoe = 0.05;     % 质量变化率权重
% gamma_qoe = 0.3;    % 初始等待时间权重
% delta_qoe = 0.2;    % 视频预期完整性权重
% 
% % 验证权重和为1
% if abs(alpha_qoe + beta_qoe + gamma_qoe + delta_qoe - 1) > 1e-6
%     error('QoE权重参数之和必须为1');
% end
% 
% %% 缓存容量
% D_eg = 30 * nf * v_chi_star; % SBS缓存容量
% D_iu = 10 * nf * v_chi_star; % IU缓存容量
% 
% %% 转码和缓冲参数
% p_sbs = 150;     % SBS计算资源，单位: Mbit/s
% p_iu = 10;      % IU计算资源，单位: Mbit/s
% D_bf = 10;     % 播放缓冲区最大容量，单位: Mbit
% t_cloud = 1;      % 云到SBS回程时间，单位: s
% t_propagation = 0.01; % SBS间传播时间，单位: s
% 
% %% 选择要优化的社区和时隙
% target_community = 1;        % 选择第1个社区进行优化
% 
% fprintf('优化目标: 社区%d\n', target_community);
% 
% %% 算法参数
% eta = 0.005;         % 步长
% epsilon_conv = 0.1; % 收敛阈值
% max_iterations = 1000; %最大次数
% run_num = 10;  % 求QoE平均值循环次数  

%% 代码运行
run("main.m");
run("Everycom_Per_IU.m");
run("IU_CacheSize.m");
run("IU_Computing.m");
run("SBS_CacheSize.m");
run("Zipf.m");


run("SBS_Computing.m");
run("Com_Parameter.m");
run("Users_Percom.m");

% run("SBS_IU_Computing.m");
% run("SBS_IU_CacheSize.m");
toc;