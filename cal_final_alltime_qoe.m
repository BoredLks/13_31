function final_alltime_qoe = cal_final_alltime_qoe(iu_count_per_community, nf, chi, Delta_t, alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, p_sbs, p_iu, D_bf, eta, epsilon_conv, max_iterations, iu_flags, user_file_req_prob, cache_decision, community_users, requested_videos, download_rates, task_assignment, initial_wait_times, r_decision, r_previous, lambda, mu_sbs, mu_iu, buffer_state)
% %% 迭代优化
% final_alltime_qoe = 0; % 初始化
% % 时间循环
% for t = 1:nf
%     convergence_lagrange = [];
%     convergence_QoE = [];
%     if t == 1
%         for i = 1:length(community_users)
%             r_previous(i, t) = chi(randi(length(chi))); % 假设视频上个视频分辨率随机
%         end
%     else
%         r_previous(:, t) = r_decision(:, t-1);
%         initial_wait_times(:)=0; % 其他时刻初始等待时间为0
%         for i = 1:length(community_users)
%             buffer_state(i, t) = max(0, min(D_bf , buffer_state(i, t-1) + (download_rates(i,t-1)-r_decision(i,t-1)) * Delta_t));
%         end
%     end
%     %% 主优化循环
%     for iter = 1:max_iterations
%         % 初始化
%         sbs_load = 0;
%         iu_load = zeros(1,iu_count_per_community);
%         lagrange = 0;
%         qoe_all_user = 0;
% 
%         for n = 1:length(community_users)
%             % 只有非最高分辨率的用户才计入负载
%             if task_assignment(n,t) == 0 && r_decision(n,t) < chi(end)
%                 sbs_load = sbs_load + r_decision(n,t);
%             end
%         end
% 
%         for n = 1:iu_count_per_community
%             for m = 1:length(community_users)
%                 % 只有非最高分辨率的用户才计入负载
%                 if task_assignment(m,t) == n && r_decision(m,t) < chi(end)
%                     iu_load(n) = iu_load(n) + r_decision(m,t);
%                 end
%             end
%         end
% 
% 
%         % 计算目标函数值
%         for i = 1:length(community_users) %对用户目标社区的第i个用户进行优化
% 
%             % 确定对应的计算资源约束的μ值
%             if task_assignment(i,t) == 0
%                 mu_val = mu_sbs(t);
%                 fuzai = sbs_load;
%                 rongliang = p_sbs;
% 
%             else
%                 mu_val = mu_iu(task_assignment(i,t),t);
%                 fuzai = iu_load(task_assignment(i,t));
%                 rongliang = p_iu;
% 
%             end
% 
%             % 如果选择最高分辨率，则μ值设为0（不考虑计算资源约束）
%             if r_decision(i,t) >= chi(end) - 1e-6
%                 mu_val = 0;
%             end
% 
%             % 计算各项QoE指标
%             VQ = r_decision(i,t);
%             SW = (r_decision(i,t) - r_previous(i,t))^2;
%             T_w = initial_wait_times(i);
% 
%             % 记录当前缓冲区状态
%             Bu_current = buffer_state(i,t);
% 
%             % 计算视频剩余时间
%             T_remaining = nf * Delta_t - (t - 1) * Delta_t;
% 
%             % 计算回程下载速率
%             bk_rate = download_rates(i,t);
% 
%             % 计算预期完整性 EI
%             if bk_rate > r_decision(i,t)
%                 EI = 1;
%             else
%                 if T_remaining > 0
%                     numerator = Bu_current + T_remaining * bk_rate;
%                     denominator = r_decision(i,t) * T_remaining;
%                     if denominator > 0
%                         EI = min(1, numerator / denominator);
%                     end
%                 else
%                     EI = 1;
%                 end
%             end
% 
%             % 计算QoE
%             qoe_i = alpha_qoe * VQ - beta_qoe * SW - gamma_qoe * T_w + delta_qoe * EI ;
% 
%             % 计算惩罚
%             violate = lambda(i,t) * (r_decision(i,t) - download_rates(i,t)) + mu_val * (fuzai-rongliang);
% 
%             % 所有用户求和
%             lagrange = lagrange +  qoe_i - violate;
%             qoe_all_user = qoe_all_user + qoe_i;
%         end
% 
%         %% 检查收敛性
%         convergence_QoE(end+1) = qoe_all_user;
%         convergence_lagrange(end+1) = lagrange;
%         if iter > 1
%             change = abs(convergence_lagrange(end) - convergence_lagrange(end-1));
%             if change < epsilon_conv
%                 break;
%             end
%         end
% 
%         %% 步骤1: 更新分辨率决策  %% 考虑把不需要计算的计算节点不作为约束进行处理。
%         for i = 1:length(community_users)
% 
%             % 计算视频剩余时间
%             T_remaining = nf * Delta_t - (t - 1) * Delta_t;
% 
%             % 计算当前缓冲区状态
%             Bu_current = buffer_state(i,t);
% 
%             % 计算回程下载速率
%             bk_rate = download_rates(i,t);
% 
%             % 重新计算预期完整性 EI
%             if bk_rate > r_decision(i,t)
%                 EI = 1;
%             else
%                 if T_remaining > 0
%                     numerator = Bu_current + T_remaining * bk_rate;
%                     denominator = r_decision(i,t) * T_remaining;
%                     if denominator > 0
%                         EI = min(1, numerator / denominator);
%                     end
%                 else
%                     EI = 1;
%                 end
%             end
% 
%             % 确定对应的计算资源约束的μ值
%             if task_assignment(i,t) == 0
%                 mu_val = mu_sbs(t);
%             else
%                 mu_val = mu_iu(task_assignment(i,t),t);
%             end
%             % 最高分辨率不受计算资源约束 =====
%             if r_decision(i,t) >= chi(end) - 1e-6
%                 mu_val = 0;
%             end
% 
%             % 根据EI值选择更新方程
%             if EI >= 0.999  % 近似等于1
%                 % 使用线性更新方程
%                 r_new = (alpha_qoe - lambda(i,t) - mu_val) / (2 * beta_qoe) + r_previous(i,t);
%             else
%                 % 使用三次方程求解
%                 % 系数计算
%                 a = 2 * beta_qoe * T_remaining;
%                 b = -(alpha_qoe + 2 * beta_qoe * r_previous(i,t) - lambda(i,t) - mu_val) * T_remaining;
%                 c = 0;
%                 d = delta_qoe * (Bu_current + T_remaining * bk_rate);
% 
%                 % 求解三次方程 ar³ + br² + cr + d = 0
%                 coeffs = [a, b, c, d];
%                 roots_complex = roots(coeffs);
% 
%                 % 找到实数根
%                 real_roots = real(roots_complex(abs(imag(roots_complex)) < 1e-10));
% 
%                 % 选择在合理范围内的根
%                 valid_roots = real_roots(real_roots > 0 & real_roots <= max(chi));
% 
%                 if ~isempty(valid_roots)
%                     r_new = valid_roots(1);
%                 else
%                     % 如果没有有效根，使用线性近似
%                     r_new = (alpha_qoe - lambda(i,t) - mu_val) / (2 * beta_qoe) + r_previous(i,t);
%                 end
%             end
% 
%             % 约束r在有效范围内
%             r_new = min(r_new, download_rates(i,t)); % 不能超过下载速率
%             r_new = max(r_new, 0); %大于零
%             r_decision(i,t) = r_new;
%         end
% 
%         %% 步骤2: 更新对偶变量
%         % 更新λ
%         for i = 1:length(community_users)
%             lambda(i,t) = max(0, lambda(i,t) + eta * (r_decision(i,t) - download_rates(i,t)));
%         end
% 
%         % 更新μ_sbs
%         % 更新μ时排除最高分辨率用户
%         sbs_load = 0;
%         iu_load = zeros(1,iu_count_per_community);
%         for n = 1:length(community_users)
%             if task_assignment(n,t) == 0 && r_decision(n,t) < chi(end)
%                 sbs_load = sbs_load + r_decision(n,t);
%             end
%         end
%         mu_sbs(t) = max(0, mu_sbs(1,t) + eta * (sbs_load - p_sbs));
% 
%         % 更新μ_iu
%         for n = 1:iu_count_per_community
%             for m = 1:length(community_users)
%                 if task_assignment(m,t) == n && r_decision(m,t) < chi(end)
%                     iu_load(n) = iu_load(n) + r_decision(m,t);
%                 end
%             end
%             mu_iu(n,t) = max(0, mu_iu(n,t) + eta * (iu_load(n) - p_iu));
%         end
% 
%     end
% 
%     %% 分辨率恢复到离散值
%     r_final = zeros(length(community_users), 1);
% 
%     % 步骤1: 初始化分辨率选择（higher or lower）
%     for i = 1:length(community_users)
% 
%         % 确定higher和lower
%         if r_decision(i,t) <= chi(1)
%             lower = chi(1);
%             higher = chi(1);
%         elseif r_decision(i,t) >= chi(end)
%             lower = chi(end);
%             higher = chi(end);
%         else
%             % 找到r_decision(i)的上下界
%             upper_idx = find(chi >= r_decision(i,t), 1, 'first');
%             if upper_idx == 1
%                 lower = chi(1);
%                 higher = chi(1);
%             else
%                 lower = chi(upper_idx - 1);
%                 higher = chi(upper_idx);
%             end
%         end
% 
%         % 计算选择lower时的QoE
%         qoe_lower = calculate_qoe(Bu_current, lower, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
% 
%         % 计算选择higher时的QoE
%         qoe_higher = calculate_qoe(Bu_current, higher, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
% 
%         % 选择QoE更高的分辨率作为初始值
%         if qoe_higher >= qoe_lower
%             r_final(i) = higher;
%         else
%             r_final(i) = lower;
%         end
%     end
% 
% 
%     % 步骤2: 检查计算资源约束并进行降级
%     % 按不同计算节点分别检查约束
% 
%     % 检查SBS约束
%     sbs_users = find(task_assignment(:,t) == 0);
%     sbs_users_need_transcode = [];
%     for idx = 1:length(sbs_users)
%         if r_final(sbs_users(idx)) < chi(end)
%             sbs_users_need_transcode(end+1) = sbs_users(idx);
%         end
%     end
%     sbs_total_load = sum(r_final(sbs_users_need_transcode));
% 
%     if sbs_total_load > p_sbs
% 
%         % 计算每个SBS用户降级的QoE损失
%         qoe_losses = [];
%         user_resolution_info = [];
% 
%         for idx = 1:length(sbs_users_need_transcode)
%             i = sbs_users_need_transcode(idx);
%             current_res = r_final(i);
% 
%             % 找到lower分辨率
%             current_res_idx = find(chi == current_res);
%             if current_res_idx > 1
%                 lower_res = chi(current_res_idx - 1);
% 
%                 % 计算当前QoE和降级后QoE
%                 qoe_current = calculate_qoe(Bu_current, current_res, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
%                 qoe_lower = calculate_qoe(Bu_current, lower_res, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
% 
%                 % 计算QoE损失
%                 qoe_loss = qoe_lower - qoe_current;
%                 qoe_losses(end+1) = qoe_loss;
%                 user_resolution_info(end+1, :) = [i, current_res, lower_res, qoe_loss];
%             end
%         end
% 
%         % 按QoE损失升序排列（优先降级损失最小的）
%         if ~isempty(qoe_losses)
%             [~, sort_indices] = sort(qoe_losses, 'descend'); % 升序排列
%             sorted_user_info = user_resolution_info(sort_indices, :);
% 
%             % 逐一降级直到满足约束
%             for idx = 1:size(sorted_user_info, 1)
%                 user_i = sorted_user_info(idx, 1);
%                 lower_res = sorted_user_info(idx, 3);
% 
%                 % 尝试降级
%                 r_final(user_i) = lower_res;
% 
%                 % 重新计算SBS负载
%                 new_sbs_load = 0;
%                 for check_idx = 1:length(sbs_users_need_transcode)
%                     check_user = sbs_users_need_transcode(check_idx);
%                     if r_final(check_user) < chi(end)
%                         new_sbs_load = new_sbs_load + r_final(check_user);
%                     end
%                 end
% 
%                 if new_sbs_load <= p_sbs
%                     break; % 满足约束，停止降级
%                 end
%             end
%         end
%     end
% 
%     % 检查每个IU的约束
%     for j = 1:iu_count_per_community
%         iu_users = find(task_assignment(:,t) == j);
%         iu_users_need_transcode = [];
%         for idx = 1:length(iu_users)
%             if r_final(iu_users(idx)) < chi(end)
%                 iu_users_need_transcode(end+1) = iu_users(idx);
%             end
%         end
% 
%         if ~isempty(iu_users_need_transcode)
%             iu_total_load = sum(r_final(iu_users_need_transcode));
% 
%             if iu_total_load > p_iu
% 
%                 % 计算每个IU用户降级的QoE损失
%                 qoe_losses = [];
%                 user_resolution_info = [];
% 
%                 for idx = 1:length(iu_users_need_transcode)
%                     i = iu_users_need_transcode(idx);
%                     current_res = r_final(i);
% 
%                     % 找到lower分辨率
%                     current_res_idx = find(chi == current_res);
%                     if current_res_idx > 1
%                         lower_res = chi(current_res_idx - 1);
% 
%                         % 计算当前QoE和降级后QoE
%                         qoe_current = calculate_qoe(Bu_current, current_res, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
%                         qoe_lower = calculate_qoe(Bu_current, lower_res, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
% 
%                         % 计算QoE损失
%                         qoe_loss = qoe_lower - qoe_current;
%                         qoe_losses(end+1) = qoe_loss;
%                         user_resolution_info(end+1, :) = [i, current_res, lower_res, qoe_loss];
%                     end
%                 end
% 
%                 % 按QoE损失升序排列（优先降级损失最小的）
%                 if ~isempty(qoe_losses)
%                     [~, sort_indices] = sort(qoe_losses, 'descend'); % 升序排列
%                     sorted_user_info = user_resolution_info(sort_indices, :);
% 
%                     % 逐一降级直到满足约束
%                     for idx = 1:size(sorted_user_info, 1)
%                         user_i = sorted_user_info(idx, 1);
%                         lower_res = sorted_user_info(idx, 3);
% 
%                         % 尝试降级
%                         r_final(user_i) = lower_res;
% 
%                         % 重新计算IU负载
%                         new_iu_load = 0;
%                         for check_idx = 1:length(iu_users_need_transcode)
%                             check_user = iu_users_need_transcode(check_idx);
%                             if r_final(check_user) < chi(end)
%                                 new_iu_load = new_iu_load + r_final(check_user);
%                             end
%                         end
% 
%                         if new_iu_load <= p_iu
%                             break; % 满足约束，停止降级
%                         end
%                     end
%                 end
%             end
%         end
%     end
% 
%     r_decision(:,t) = r_final; % 赋值为最终确定的在t时刻时的r_decision
% 
%     %% 计算t时刻时的最终QoE
%     final_total_qoe = 0;
%     for i = 1:length(community_users)
%         % user_idx = community_users(i);
%         % requested_video = requested_videos(i);
% 
%         % 计算最终QoE
%         qoe_i = calculate_qoe(Bu_current, r_final(i), r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
%         final_total_qoe = final_total_qoe + qoe_i;  %总QoE值
%         % final_total_qoe = final_total_qoe + user_file_req_prob(user_idx, requested_video) * qoe_i;  %加权总QoE值
% 
%     end
% 
%     final_alltime_qoe = final_alltime_qoe + final_total_qoe;
% end





























%% 迭代优化
final_alltime_qoe = 0; % 初始化
% 时间循环
for t = 1:nf
    convergence_lagrange = [];
    convergence_QoE = [];
    if t == 1
        for i = 1:length(community_users)
            r_previous(i, t) = chi(randi(length(chi))); % 假设视频上个视频分辨率随机
        end
    else
        r_previous(:, t) = r_decision(:, t-1);
        initial_wait_times(:)=0; % 其他时刻初始等待时间为0
        for i = 1:length(community_users)
            buffer_state(i, t) = max(0, min(D_bf , buffer_state(i, t-1) + (download_rates(i,t-1)-r_decision(i,t-1)) * Delta_t));
        end
    end
    %% 主优化循环
    for iter = 1:max_iterations
        % 初始化
        sbs_load = 0;
        iu_load = zeros(1,iu_count_per_community);
        lagrange = 0;
        qoe_all_user = 0;

        for n = 1:length(community_users)
            % 只有非最高分辨率的用户才计入负载
            if task_assignment(n,t) == 0 && r_decision(n,t) < chi(end)
                sbs_load = sbs_load + r_decision(n,t);
            end
        end

        for n = 1:iu_count_per_community
            for m = 1:length(community_users)
                % 只有非最高分辨率的用户才计入负载
                if task_assignment(m,t) == n && r_decision(m,t) < chi(end)
                    iu_load(n) = iu_load(n) + r_decision(m,t);
                end
            end
        end


        % 计算目标函数值
        for i = 1:length(community_users) %对用户目标社区的第i个用户进行优化

            % 确定对应的计算资源约束的μ值
            if task_assignment(i,t) == 0
                mu_val = mu_sbs(t);
                fuzai = sbs_load;
                rongliang = p_sbs;

            else
                mu_val = mu_iu(task_assignment(i,t),t);
                fuzai = iu_load(task_assignment(i,t));
                rongliang = p_iu;

            end

            % 如果选择最高分辨率，则μ值设为0（不考虑计算资源约束）
            if r_decision(i,t) >= chi(end) - 1e-6
                mu_val = 0;
            end

            % 计算各项QoE指标
            VQ = r_decision(i,t);
            SW = (r_decision(i,t) - r_previous(i,t))^2;
            T_w = initial_wait_times(i);

            % 记录当前缓冲区状态
            Bu_current = buffer_state(i,t);

            % 计算视频剩余时间
            T_remaining = nf * Delta_t - (t - 1) * Delta_t;

            % 计算回程下载速率
            bk_rate = download_rates(i,t);

            % 计算预期完整性 EI
            if bk_rate > r_decision(i,t)
                EI = 1;
            else
                if T_remaining > 0
                    numerator = Bu_current + T_remaining * bk_rate;
                    denominator = r_decision(i,t) * T_remaining;
                    if denominator > 0
                        EI = min(1, numerator / denominator);
                    end
                else
                    EI = 1;
                end
            end

            % 计算QoE
            qoe_i = alpha_qoe * VQ - beta_qoe * SW - gamma_qoe * T_w + delta_qoe * EI ;

            % 计算惩罚
            violate = lambda(i,t) * (r_decision(i,t) - download_rates(i,t)) + mu_val * (fuzai-rongliang);

            % 所有用户求和
            lagrange = lagrange +  qoe_i - violate;
            qoe_all_user = qoe_all_user + qoe_i;
        end

        %% 检查收敛性
        convergence_QoE(end+1) = qoe_all_user;
        convergence_lagrange(end+1) = lagrange;
        if iter > 1
            change = abs(convergence_lagrange(end) - convergence_lagrange(end-1));
            if change < epsilon_conv
                break;
            end
        end

        %% 步骤1: 更新分辨率决策  考虑把不需要计算的计算节点不作为约束进行处理。
        for i = 1:length(community_users)
            
            user_idx = community_users(i);
            if iu_flags(user_idx) == 1 && cache_decision(user_idx, requested_videos(i)) == 1
                % IU用户缓存命中，直接使用最高分辨率并跳过此用户的更新
                r_decision(i, t) = chi(end);
                continue;
            end

            % 计算视频剩余时间
            T_remaining = nf * Delta_t - (t - 1) * Delta_t;

            % 计算当前缓冲区状态
            Bu_current = buffer_state(i,t);

            % 计算回程下载速率
            bk_rate = download_rates(i,t);

            % 重新计算预期完整性 EI
            if bk_rate > r_decision(i,t)
                EI = 1;
            else
                if T_remaining > 0
                    numerator = Bu_current + T_remaining * bk_rate;
                    denominator = r_decision(i,t) * T_remaining;
                    if denominator > 0
                        EI = min(1, numerator / denominator);
                    end
                else
                    EI = 1;
                end
            end

            % 确定对应的计算资源约束的μ值
            if task_assignment(i,t) == 0
                mu_val = mu_sbs(t);
            else
                mu_val = mu_iu(task_assignment(i,t),t);
            end
            % 最高分辨率不受计算资源约束 =====
            if r_decision(i,t) >= chi(end) - 1e-6
                mu_val = 0;
            end

            % 根据EI值选择更新方程
            if EI >= 0.999  % 近似等于1
                % 使用线性更新方程
                r_new = (alpha_qoe - lambda(i,t) - mu_val) / (2 * beta_qoe) + r_previous(i,t);
            else
                % 使用三次方程求解
                % 系数计算
                a = 2 * beta_qoe * T_remaining;
                b = -(alpha_qoe + 2 * beta_qoe * r_previous(i,t) - lambda(i,t) - mu_val) * T_remaining;
                c = 0;
                d = delta_qoe * (Bu_current + T_remaining * bk_rate);

                % 求解三次方程 ar³ + br² + cr + d = 0
                coeffs = [a, b, c, d];
                roots_complex = roots(coeffs);

                % 找到实数根
                real_roots = real(roots_complex(abs(imag(roots_complex)) < 1e-10));

                % 选择在合理范围内的根
                valid_roots = real_roots(real_roots > 0 & real_roots <= max(chi));

                if ~isempty(valid_roots)
                    r_new = valid_roots(1);
                else
                    % 如果没有有效根，使用线性近似
                    r_new = (alpha_qoe - lambda(i,t) - mu_val) / (2 * beta_qoe) + r_previous(i,t);
                end
            end

            % 约束r在有效范围内
            r_new = min(r_new, download_rates(i,t)); % 不能超过下载速率
            r_new = max(r_new, 0); %大于零
            r_decision(i,t) = r_new;
        end

        %% 步骤2: 更新对偶变量
        % 更新λ
        for i = 1:length(community_users)
            lambda(i,t) = max(0, lambda(i,t) + eta * (r_decision(i,t) - download_rates(i,t)));
        end

        % 更新μ_sbs
        % 更新μ时排除最高分辨率用户
        sbs_load = 0;
        iu_load = zeros(1,iu_count_per_community);
        for n = 1:length(community_users)
            if task_assignment(n,t) == 0 && r_decision(n,t) < chi(end)
                sbs_load = sbs_load + r_decision(n,t);
            end
        end
        mu_sbs(t) = max(0, mu_sbs(1,t) + eta * (sbs_load - p_sbs));

        % 更新μ_iu
        for n = 1:iu_count_per_community
            for m = 1:length(community_users)
                if task_assignment(m,t) == n && r_decision(m,t) < chi(end)
                    iu_load(n) = iu_load(n) + r_decision(m,t);
                end
            end
            mu_iu(n,t) = max(0, mu_iu(n,t) + eta * (iu_load(n) - p_iu));
        end

    end

    %% 分辨率恢复到离散值
    r_final = zeros(length(community_users), 1);

    % 步骤1: 初始化分辨率选择（higher or lower）
    for i = 1:length(community_users)

        % 确定higher和lower
        if r_decision(i,t) <= chi(1)
            lower = chi(1);
            higher = chi(1);
        elseif r_decision(i,t) >= chi(end)
            lower = chi(end);
            higher = chi(end);
        else
            % 找到r_decision(i)的上下界
            upper_idx = find(chi >= r_decision(i,t), 1, 'first');
            if upper_idx == 1
                lower = chi(1);
                higher = chi(1);
            else
                lower = chi(upper_idx - 1);
                higher = chi(upper_idx);
            end
        end

        % 计算选择lower时的QoE
        qoe_lower = calculate_qoe(Bu_current, lower, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);

        % 计算选择higher时的QoE
        qoe_higher = calculate_qoe(Bu_current, higher, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);

        % 选择QoE更高的分辨率作为初始值
        if qoe_higher >= qoe_lower
            r_final(i) = higher;
        else
            r_final(i) = lower;
        end
    end


    % 步骤2: 检查计算资源约束并进行降级
    % 按不同计算节点分别检查约束

    % 检查SBS约束
    sbs_users = find(task_assignment(:,t) == 0);
    sbs_users_need_transcode = [];
    for idx = 1:length(sbs_users)
        if r_final(sbs_users(idx)) < chi(end)
            sbs_users_need_transcode(end+1) = sbs_users(idx);
        end
    end
    sbs_total_load = sum(r_final(sbs_users_need_transcode));

    if sbs_total_load > p_sbs

        % 计算每个SBS用户降级的QoE损失
        qoe_losses = [];
        user_resolution_info = [];

        for idx = 1:length(sbs_users_need_transcode)
            i = sbs_users_need_transcode(idx);
            current_res = r_final(i);

            % 找到lower分辨率
            current_res_idx = find(chi == current_res);
            if current_res_idx > 1
                lower_res = chi(current_res_idx - 1);

                % 计算当前QoE和降级后QoE
                qoe_current = calculate_qoe(Bu_current, current_res, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
                qoe_lower = calculate_qoe(Bu_current, lower_res, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);

                % 计算QoE损失
                qoe_loss = qoe_lower - qoe_current;
                qoe_losses(end+1) = qoe_loss;
                user_resolution_info(end+1, :) = [i, current_res, lower_res, qoe_loss];
            end
        end

        % 按QoE损失升序排列（优先降级损失最小的）
        if ~isempty(qoe_losses)
            [~, sort_indices] = sort(qoe_losses, 'descend'); % 升序排列
            sorted_user_info = user_resolution_info(sort_indices, :);

            % 逐一降级直到满足约束
            for idx = 1:size(sorted_user_info, 1)
                user_i = sorted_user_info(idx, 1);
                lower_res = sorted_user_info(idx, 3);

                % 尝试降级
                r_final(user_i) = lower_res;

                % 重新计算SBS负载
                new_sbs_load = 0;
                for check_idx = 1:length(sbs_users_need_transcode)
                    check_user = sbs_users_need_transcode(check_idx);
                    if r_final(check_user) < chi(end)
                        new_sbs_load = new_sbs_load + r_final(check_user);
                    end
                end

                if new_sbs_load <= p_sbs
                    break; % 满足约束，停止降级
                end
            end
        end
    end

    % 检查每个IU的约束
    for j = 1:iu_count_per_community
        iu_users = find(task_assignment(:,t) == j);
        iu_users_need_transcode = [];
        for idx = 1:length(iu_users)
            if r_final(iu_users(idx)) < chi(end)
                iu_users_need_transcode(end+1) = iu_users(idx);
            end
        end

        if ~isempty(iu_users_need_transcode)
            iu_total_load = sum(r_final(iu_users_need_transcode));

            if iu_total_load > p_iu

                % 计算每个IU用户降级的QoE损失
                qoe_losses = [];
                user_resolution_info = [];

                for idx = 1:length(iu_users_need_transcode)
                    i = iu_users_need_transcode(idx);
                    current_res = r_final(i);

                    % 找到lower分辨率
                    current_res_idx = find(chi == current_res);
                    if current_res_idx > 1
                        lower_res = chi(current_res_idx - 1);

                        % 计算当前QoE和降级后QoE
                        qoe_current = calculate_qoe(Bu_current, current_res, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
                        qoe_lower = calculate_qoe(Bu_current, lower_res, r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);

                        % 计算QoE损失
                        qoe_loss = qoe_lower - qoe_current;
                        qoe_losses(end+1) = qoe_loss;
                        user_resolution_info(end+1, :) = [i, current_res, lower_res, qoe_loss];
                    end
                end

                % 按QoE损失升序排列（优先降级损失最小的）
                if ~isempty(qoe_losses)
                    [~, sort_indices] = sort(qoe_losses, 'descend'); % 升序排列
                    sorted_user_info = user_resolution_info(sort_indices, :);

                    % 逐一降级直到满足约束
                    for idx = 1:size(sorted_user_info, 1)
                        user_i = sorted_user_info(idx, 1);
                        lower_res = sorted_user_info(idx, 3);

                        % 尝试降级
                        r_final(user_i) = lower_res;

                        % 重新计算IU负载
                        new_iu_load = 0;
                        for check_idx = 1:length(iu_users_need_transcode)
                            check_user = iu_users_need_transcode(check_idx);
                            if r_final(check_user) < chi(end)
                                new_iu_load = new_iu_load + r_final(check_user);
                            end
                        end
                        
                        if new_iu_load > p_iu
                            r_final(user_i) = chi(find(chi == lower_res) - 1); % 如果还是不满足约束，IU的连接用户还可以再降一级
                        else
                            break; % 满足约束，停止降级
                        end
                    end
                end
            end
        end
    end

    r_decision(:,t) = r_final; % 赋值为最终确定的在t时刻时的r_decision

    %% 计算t时刻时的最终QoE
    final_total_qoe = 0;
    for i = 1:length(community_users)
        % user_idx = community_users(i);
        % requested_video = requested_videos(i);

        % 计算最终QoE
        qoe_i = calculate_qoe(Bu_current, r_final(i), r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
        final_total_qoe = final_total_qoe + qoe_i;  %总QoE值
        % final_total_qoe = final_total_qoe + user_file_req_prob(user_idx, requested_video) * qoe_i;  %加权总QoE值

    end
    
    final_alltime_qoe = final_alltime_qoe + final_total_qoe;
end













end
