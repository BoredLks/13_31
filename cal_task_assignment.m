function task_assignment = cal_task_assignment(iu_count_per_community, nf, iu_coverage, target_community, user_positions, iu_indices, cache_decision, community_users, requested_videos)
%% 计算计算卸载分配
task_assignment = zeros(length(community_users), nf); % 初始化

for t = 1:nf
    for i = 1:length(community_users)
        user_idx = community_users(i);
        requested_video = requested_videos(i);

        assigned = false;
        for j = 1:iu_count_per_community
            iu_idx = iu_indices(target_community, j);
            if cache_decision(iu_idx, requested_video) == 1
                dist = sqrt(sum((user_positions(user_idx, :, t) - user_positions(iu_idx, :, t)).^2));

                if dist <= iu_coverage
                    task_assignment(i,t) = j;
                    assigned = true;
                    break;
                end
            end
        end

        if ~assigned
            task_assignment(i,t) = 0;
        end
    end
end

end