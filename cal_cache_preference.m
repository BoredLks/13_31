function cache_preference = cal_cache_preference(h, total_users, F, user_file_req_prob, joint_edge_weights)
%% 计算缓存偏好
cache_preference = zeros(total_users + h, F);

% 计算每个节点对每个文件的缓存偏好
for i = 1:(total_users + h)
    
    out_neighbors = find(joint_edge_weights(i, :) > 0);    % 找出所有能接收i提供服务的节点（出度邻居）
    
    out_neighbors = out_neighbors(~ismember(out_neighbors, i)); %排除IU节点自己的偏好。

    if ~isempty(out_neighbors)
        for k = 1:F
            preference_sum = 0;
            
            for j = out_neighbors
                if j <= total_users % 只考虑用户节点的请求概率
                    preference_sum = preference_sum + user_file_req_prob(j, k) * joint_edge_weights(i, j);
                end
            end
            
            cache_preference(i, k) = preference_sum / length(out_neighbors);
        end
    end
end
end