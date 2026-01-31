function requested_videos = cal_requested_videos(user_file_req_prob, community_users)
%% 获取请求视频
requested_videos = zeros(length(community_users), 1);
for i = 1:length(community_users)
    user_idx = community_users(i);
    % 根据用户的文件请求概率分布随机选择一个视频
    cumulative_prob = cumsum(user_file_req_prob(user_idx, :));
    rand_val = rand();
    requested_videos(i) = find(cumulative_prob >= rand_val, 1);
end

end
