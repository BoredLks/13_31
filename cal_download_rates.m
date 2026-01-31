function download_rates = cal_download_rates(h, iu_count_per_community, nf, sbs_coverage, iu_coverage, B, P_sbs, P_iu, N0, K, epsilon, slowfading_dB, target_community, sbs_positions, user_positions, iu_indices, iu_flags, cache_decision, community_users, requested_videos, task_assignment)
%% 计算每个时隙的用户下载速率和任务分配
download_rates = zeros(length(community_users), nf); % 初始化
R_max = B * log2(1 + (P_sbs * K * 1 * 1) / N0); % 计算最大理论速率

% 对每个时隙计算下载速率和任务分配
for t = 1:nf
    for i = 1:length(community_users)
        user_idx = community_users(i);
        requested_video = requested_videos(i);

        % 判断用户下载速率的情况
        if iu_flags(user_idx) == 1 && cache_decision(user_idx, requested_video) == 1 % Case 1: IU用户自己请求且缓存了请求的文件
            download_rates(i, t) = R_max;
        else 
            % 寻找最佳的服务节点
            best_rate = 0;

            % 检查社区内其他IU是否缓存了该文件 % Case 2: 附近IU缓存了请求的文件
            for j = 1:iu_count_per_community
                iu_idx = iu_indices(target_community, j);
                if cache_decision(iu_idx, requested_video) == 1
                    % 计算距离
                    dist = sqrt(sum((user_positions(user_idx, :, t) - user_positions(iu_idx, :, t)).^2));

                    if dist <= iu_coverage
                        % 计算信道增益和干扰
                        % 大尺度
                        slowfading_sd = slowfading_dB / (10 * log10(exp(1)));% 将 dB 转换为对数后的标准差
                        slowfading_avg = -slowfading_sd^2 / 2;% 假设均值为 0 dB，则对数后的均值为 mu = -sigma^2 / 2
                        vartheta = lognrnd(slowfading_avg, slowfading_sd);
                        % 小尺度
                        re_real = randn();
                        re_imag = randn();
                        re = sqrt(re_real + 1i * re_imag);
                        xi = abs(re) ^ 2;
                        % 增益
                        channel_gain = K * vartheta * xi * dist^(-epsilon);

                        % 计算干扰
                        interference = 0;
                        for k = 1:iu_count_per_community
                            other_iu_idx = iu_indices(target_community, k);
                            if other_iu_idx ~= iu_idx && other_iu_idx ~= user_idx
                                dist_interferer = sqrt(sum((user_positions(other_iu_idx, :, t) - user_positions(user_idx, :, t)).^2));

                                if dist_interferer <= iu_coverage && ~isempty(find(task_assignment(:,t) == find(iu_indices(target_community, :) == other_iu_idx), 1))
                                    
                                    % 大尺度
                                    vartheta_int = lognrnd(slowfading_avg, slowfading_sd);
                                    % 小尺度
                                    re_real_int = randn();
                                    re_imag_int = randn();
                                    re_int = sqrt(re_real_int + 1i * re_imag_int);
                                    xi_int = abs(re_int) ^ 2;
                                    % 总干扰
                                    interference_gain = K * vartheta_int * xi_int * dist_interferer^(-epsilon);
                                    interference = interference + P_iu * interference_gain;
                                end
                            end
                        end

                        rate = B * log2(1 + (P_iu * channel_gain) / (N0 + interference));

                        best_rate = max(best_rate, rate);
                    end
                end
            end

           if (best_rate <= 0.5) % 除非IU的速率实在是太低，那么只能选择SBS进行传输
                best_rate = 0;
           end

            % 如果没有找到合适的IU，使用SBS % Case 3: 所处区域的SBS传输，缓存可能在本区域SBS，云，其他SBS上
            if best_rate == 0
                dist_to_sbs = sqrt(sum((user_positions(user_idx, :, t) - sbs_positions(target_community, :)).^2));

                if dist_to_sbs <= sbs_coverage
                    % 计算信道增益和干扰
                    % 大尺度
                    slowfading_sd = slowfading_dB / (10 * log10(exp(1)));% 将 dB 转换为对数后的标准差
                    slowfading_avg = -slowfading_sd^2 / 2;% 假设均值为 0 dB，则对数后的均值为 mu = -sigma^2 / 2
                    vartheta = lognrnd(slowfading_avg, slowfading_sd);
                    % 小尺度
                    re_real = randn();
                    re_imag = randn();
                    re = sqrt(re_real + 1i * re_imag);
                    xi = abs(re) ^ 2;
                    % 增益
                    channel_gain = K * vartheta * xi * dist_to_sbs^(-epsilon);

                    % 计算干扰
                    interference = 0;
                    for m = 1:h
                        if m ~= target_community
                             % 计算其他SBS到用户的距离
                            dist_interferer_sbs = sqrt(sum((user_positions(user_idx, :, t) - sbs_positions(m, :)).^2));
                            % 如果其他SBS在用户的接收范围内（可能造成干扰）%if dist_interferer_sbs <= sbs_coverage %end
                            % 大尺度
                            vartheta_int = lognrnd(slowfading_avg, slowfading_sd);
                            % 小尺度
                            re_real_int = randn();
                            re_imag_int = randn();
                            re_int = sqrt(re_real_int + 1i * re_imag_int);
                            xi_int = abs(re_int) ^ 2;
                            % 总干扰
                            interference_gain = K * vartheta_int * xi_int * dist_interferer_sbs^(-epsilon);
                            interference = interference + P_sbs * interference_gain;

                        end
                    end

                    best_rate = B * log2(1 + (P_sbs * channel_gain) / (N0 + interference));
                end
            end

            download_rates(i, t) = best_rate;
        end
    end
end
end