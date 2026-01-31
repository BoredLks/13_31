function [iu_indices,iu_flags] = cal_Select_iu(h, user_per_community, total_users, iu_count_per_community, iu_coverage, user_positions)

% %% 选择每个社区的重要用户作为缓存节点
% % 目标：选择的IU覆盖更多人数，且覆盖范围尽量不重叠
% iu_indices = zeros(h, iu_count_per_community);
% iu_flags = zeros(total_users, 1); % 标记哪些用户是IU
% 
% % 使用过去时间的平均位置来计算用户密度和覆盖范围
% initial_positions = mean(user_positions, 3);
% 
% for m = 1:h
%     start_idx = (m-1) * user_per_community + 1;
%     end_idx = m * user_per_community;
%     community_users = start_idx:end_idx;
% 
%     % 使用贪心算法选择IU，最大化覆盖用户数且最小化重叠
%     selected_ius = [];
%     covered_users = [];
% 
%     for iu_round = 1:iu_count_per_community
%         best_candidate = -1;
%         best_score = -inf;
%         best_new_coverage = [];
% 
%         % 遍历所有候选用户
%         for candidate_idx = 1:length(community_users)
%             candidate_user = community_users(candidate_idx);
% 
%             % 跳过已选择的IU
%             if ismember(candidate_user, selected_ius)
%                 continue;
%             end
% 
%             candidate_pos = initial_positions(candidate_user, :);
% 
%             % 计算该候选用户的覆盖范围内的用户
%             coverage_users = [];
%             for j = 1:length(community_users)
%                 other_user = community_users(j);
%                 other_pos = initial_positions(other_user, :);
% 
%                 % 计算距离
%                 distance = sqrt(sum((candidate_pos - other_pos).^2));
% 
%                 % 如果在D2D覆盖范围内，加入覆盖列表
%                 if distance <= iu_coverage
%                     coverage_users = [coverage_users, other_user];
%                 end
%             end
% 
%             % 计算新增覆盖的用户数（去除已被其他IU覆盖的用户）
%             new_coverage = setdiff(coverage_users, covered_users);
%             new_coverage_count = length(new_coverage);
% 
%             % 计算与已选IU的重叠程度
%             overlap_penalty = 0;
%             for selected_iu = selected_ius
%                 selected_pos = initial_positions(selected_iu, :);
%                 center_distance = sqrt(sum((candidate_pos - selected_pos).^2));
% 
%                 % 如果两个IU的覆盖圆心距离小于2倍覆盖半径，存在重叠
%                 if center_distance < 2 * iu_coverage
%                     % 计算重叠面积的近似惩罚
%                     if center_distance < iu_coverage
%                         % 高度重叠
%                         overlap_penalty = overlap_penalty + 1.0;
%                     else
%                         % 部分重叠，惩罚与距离成反比
%                         overlap_penalty = overlap_penalty + (2 * iu_coverage - center_distance) / iu_coverage;
%                     end
%                 end
%             end
% 
%             % 计算综合评分：新增覆盖用户数 - 重叠惩罚
%             score = new_coverage_count - 0.5 * overlap_penalty;
% 
%             % 选择最佳候选
%             if score > best_score
%                 best_score = score;
%                 best_candidate = candidate_user;
%                 best_new_coverage = new_coverage;
%             end
%         end
% 
%         % 添加最佳候选到IU列表
%         if best_candidate ~= -1
%             selected_ius = [selected_ius, best_candidate];
%             covered_users = union(covered_users, best_new_coverage);
%         else
%             break;
%         end
%     end
% 
%     % 如果选择的IU数量不足，用剩余用户按密度填充
%     while length(selected_ius) < iu_count_per_community
%         remaining_users = setdiff(community_users, selected_ius);
%         if isempty(remaining_users)
%             break;
%         end
% 
%         % 从剩余用户中选择密度最高的
%         best_density = -1;
%         best_user = -1;
% 
%         for user_idx = 1:length(remaining_users)
%             user = remaining_users(user_idx);
%             user_pos = initial_positions(user, :);
% 
%             % 计算邻居密度
%             neighbor_count = 0;
%             for j = 1:length(community_users)
%                 other_user = community_users(j);
%                 if user ~= other_user
%                     other_pos = initial_positions(other_user, :);
%                     distance = sqrt(sum((user_pos - other_pos).^2));
%                     if distance <= iu_coverage
%                         neighbor_count = neighbor_count + 1;
%                     end
%                 end
%             end
% 
%             if neighbor_count > best_density
%                 best_density = neighbor_count;
%                 best_user = user;
%             end
%         end
% 
%         if best_user ~= -1
%             selected_ius = [selected_ius, best_user];
%         else
%             break;
%         end
%     end
% 
%     % 存储选择的IU索引
%     iu_indices(m, 1:length(selected_ius)) = selected_ius;
% 
%     % 标记这些用户为IU
%     iu_flags(selected_ius) = 1;
% end







%% 随机选择每个社区的重要用户作为缓存节点
iu_indices = zeros(h, iu_count_per_community);
iu_flags = zeros(total_users, 1); % 标记哪些用户是IU
for m = 1:h
    community_user_indices = (m-1)*user_per_community + (1:user_per_community);
    selected_iu = randsample(community_user_indices, iu_count_per_community);
    iu_indices(m,:) = selected_iu;
    iu_flags(selected_iu) = 1;
end








end