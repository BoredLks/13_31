function [final_alltime_qoe,cache_hit_rate,mean_wait_time] = algorithm_5(h, user_per_community, total_users, F, iu_count_per_community, nf, chi, Delta_t, v_chi_star, region_size, sbs_coverage, iu_coverage, community_radius, max_movement_dist, D_eg, T_small, gamma_m, alpha, beta, gamma, B, P_sbs, P_iu, N0, K, epsilon, slowfading_dB, alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, p_sbs, p_iu, D_bf, t_cloud, t_propagation, target_community, eta, epsilon_conv, max_iterations, data)
%% 算法五：流行度非协助算法

[sbs_positions,per_user_positions,now_user_positions,user_community] = cal_Initial_position(h, user_per_community, total_users, region_size, community_radius, max_movement_dist, T_small); % 初始化社区位置
[iu_indices,iu_flags] = cal_Select_iu(h, user_per_community, total_users, iu_count_per_community, iu_coverage, per_user_positions); % 选择IU
[sim_matrix,user_file_req_prob] = cal_Movielens(total_users, F, gamma_m, data); % Movielens计算
sg_edge_weights = cal_sg_edge_weights(h, total_users, sbs_coverage, alpha, beta, gamma, sbs_positions, per_user_positions, user_community, sim_matrix); % 构建社交感知图
% pl_edge_weights = cal_pl_edge_weights(h, total_users, sbs_coverage, iu_coverage, T_small, B, P_sbs, P_iu, N0, K, epsilon, slowfading_dB, sbs_positions, per_user_positions, iu_flags); % 构建物理连接图
pl_edge_weights = cal_pl_edge_weights(h, total_users, iu_count_per_community, sbs_coverage, iu_coverage, T_small, B, P_sbs, P_iu, N0, K, epsilon, slowfading_dB, sbs_positions, per_user_positions, iu_indices, iu_flags);
joint_edge_weights = pl_edge_weights.*sg_edge_weights;% 构建联合图
cache_preference = cal_cache_preference(h, total_users, F, user_file_req_prob, joint_edge_weights); % 缓存偏好计算
cache_decision = cal_cache_decision_5(h, total_users, F, nf, v_chi_star, D_eg, cache_preference); % 缓存优化算法

community_users = cal_target_community_users(user_per_community, target_community); % 选择目标用户
requested_videos = cal_requested_videos(user_file_req_prob, community_users); % 请求生成
[download_rates,task_assignment] = cal_download_rates_task_assignment(h, iu_count_per_community, nf, sbs_coverage, iu_coverage, B, P_sbs, P_iu, N0, K, epsilon, slowfading_dB, target_community, sbs_positions, now_user_positions, iu_indices, iu_flags, cache_decision, community_users, requested_videos); % 对每个时隙计算下载速率和任务分配
initial_wait_times = cal_initial_wait_times(h, total_users, iu_count_per_community, v_chi_star, iu_coverage, p_sbs, p_iu, D_bf, t_cloud, t_propagation, target_community, now_user_positions, iu_indices, iu_flags, cache_decision, community_users, requested_videos, download_rates); % 计算初始等待时间（可以先计算）
mean_wait_time = mean(initial_wait_times);
[r_decision,r_previous,lambda,mu_sbs,mu_iu,buffer_state] = cal_Initialize_optimization_variables(chi, iu_count_per_community, nf, D_bf, iu_flags, cache_decision, community_users, requested_videos); % 初始化kkt参数
final_alltime_qoe = cal_final_alltime_qoe(iu_count_per_community, nf, chi, Delta_t, alpha_qoe, beta_qoe, gamma_qoe, delta_qoe, p_sbs, p_iu, D_bf, eta, epsilon_conv, max_iterations, iu_flags, user_file_req_prob, cache_decision, community_users, requested_videos, download_rates, task_assignment, initial_wait_times, r_decision, r_previous, lambda, mu_sbs, mu_iu, buffer_state); % 计算最终QoE
cache_hit_rate = cal_hit_rate(total_users, nf, target_community, cache_decision, community_users, requested_videos, task_assignment);
end