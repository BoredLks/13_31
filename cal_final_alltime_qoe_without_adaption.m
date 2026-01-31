function final_alltime_qoe = cal_final_alltime_qoe_without_adaption(nf, chi, Delta_t, alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, D_bf, community_users, download_rates, initial_wait_times, r_decision, r_previous, buffer_state)
final_alltime_qoe = 0; % 初始化总QoE值
r_final = zeros(length(community_users), 1); % 初始化最终分辨率选择
%% 时间循环
for t = 1:nf
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

    % 记录当前缓冲区状态
    Bu_current = buffer_state(i,t);

    r_final(:) = chi(2);
    r_decision(:,t) = r_final; % 赋值为最终确定的在t时刻时的r_decision

    %% 计算t时刻时的最终QoE
    final_total_qoe = 0;
    for i = 1:length(community_users)

        % 计算最终QoE
        qoe_i = calculate_qoe(Bu_current, r_final(i), r_previous(i,t), initial_wait_times(i), download_rates(i,t), alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, nf, Delta_t, t);
        final_total_qoe = final_total_qoe + qoe_i;  %总QoE值

    end
    
    final_alltime_qoe = final_alltime_qoe + final_total_qoe; %总QoE值
end

end
