function cache_hit_rate = cal_hit_rate(total_users, nf, target_community, cache_decision, community_users, requested_videos, task_assignment)

%% 计算缓存命中率
cache_hit_allcount = zeros(1, nf);
total_requests = length(requested_videos);


for t=1:nf
    for i = 1:length(community_users)
        requested_video = requested_videos(i);

        % 检查是否命中IU缓存
        if task_assignment(i, t) ~= 0
            cache_hit_allcount(t) = cache_hit_allcount(t) + 1;
            continue;
        end

        % 检查是否命中SBS缓存
        if cache_decision(total_users + target_community, requested_video) == 1
            cache_hit_allcount(t) = cache_hit_allcount(t) + 1;
        end

    end
end

cache_hit_rate = (sum(cache_hit_allcount) / total_requests) / nf;


end