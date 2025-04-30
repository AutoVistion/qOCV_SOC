%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generate SOC - qOCV Curve from Raw Data.
% 1. Find discharging event within 0.05C with constant current phase duration
% over 5 seconds.
% 2. Calculate R from each Discharge Event and save as struct
% 3. Generate SOC qOCV Curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; 

%% Directory
clc; clear; close all;

dataDir  = 'D:\JCW\KENTECH\Projects\KEPCO\ESS_Data_Preprocessing';
yearList = {'2021'}; %, '2022', '2023'};
saveDir  = fullfile(dataDir, 'qOCV_SOC\ver01');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameter Setting
Cnom                     = 1024;      % [Ah]
C_rate_limit             = Cnom*0.1;  % 0.1C Limit for Discharge Events
idle_threshold           = 1; %Cnom*0.001;

min_seg_length           = 3;         % C-rate limit 최소 구간 길이(초)
min_dSOC                 = 1;         % mininum SOC for meaningful Voltage variation
SOC_bins                 = 0:1:100;   % SOC binning

% min_thruput              = 0.05;      % charge/discharge thruput filter for SOC variation
% qOCV_limit               = Cnom*0.05;
% cc_variation_limit       = 0.05;      % for constant current detection(2%)
% min_deltaI               = 3;         % 0.05*Cnom; % minimum C-rate change

% qOCV_binned = struct();
% DCR_lookup = struct();

Vcorrected = struct();

%% Yearly SOC bin
for y = 1:length(yearList)
    yName = matlab.lang.makeValidName(yearList{y});
    for b = 1:length(SOC_bins)-1
        bName = ['bin' num2str(SOC_bins(b))]; 
        qOCV_binned.(yName).(bName) = [];
    end
end

% total_segments = 0; 
% total_steps = 0; 
% cc_filter_count = 0; 
% dynamic_filter_count = 0; 
% thruput_filter_count = 0; 
% deltaI_filter_count = 0; 
% valid_events = 0;

%% Folder Traversal
for y = 1:length(yearList)
    year = yearList{y};
    yearPath = fullfile(dataDir, year);
    monthDirs = dir(fullfile(yearPath, '20*'));

    for m = 1:length(monthDirs)
        if ~monthDirs(m).isdir, continue; end
        monthPath = fullfile(yearPath, monthDirs(m).name);
        matFiles = dir(fullfile(monthPath, 'Raw_*.mat'));

        for f = 1:length(matFiles)
            matFilePath = fullfile(monthPath, matFiles(f).name);
            [~, name, ~] = fileparts(matFiles(f).name);
            dateStr = extractAfter(name, 'Raw_');
            safeDateStr = matlab.lang.makeValidName(dateStr);
            dayField = ['dischargeEvent_', dateStr];

            % total_segments = 0;
            % total_steps = 0;
            % cc_filter_count = 0;
            % dynamic_filter_count = 0;
            % thruput_filter_count = 0;
            % deltaI_filter_count = 0;
            % valid_events = 0;

            load(matFilePath);
            t = Raw.sync_Time;
            I = Raw.Online_DC_Current;
            V = Raw.Total_Average_CV_Sum;
            soc = Raw.Total_Average_SOC;
            T_batt = Raw.Plc_Battery_Temperature_sync;
            
            eventCount = 0;
            i = 2;
            while i <= length(I) - min_seg_length
                % Discharge event starting condition: i-1: idle, i: < 0                
                if I(i) < 0 && abs(I(i-1)) <= idle_threshold
                    fprintf('[%s] 방전이벤트 후보: i = %d, I = %.2f, I(i-1) = %.2f\n', dateStr, i, I(i), I(i-1));
                    k = i - 1;
                    while k >= 1 && abs(I(k)) <= idle_threshold
                        k = k - 1;
                    end
                    idle_start = k + 1;
                    idle_range = idle_start:i-1;

                    exact_zero_idx = find(I(idle_range) == 0, 1, 'last');

                    % Set discharge start idx
                    % idle 구간에서 마지막 0
                    if ~isempty(exact_zero_idx)
                        start_idx = idle_range(exact_zero_idx);
                    else
                        start_idx = idle_range(end);
                    end

                    % Discharge end point detection
                    j = i+1;
                    while j <= length(I)
                        if I(j) >= 0 || I(j) > I(j-1)
                            break;
                        end
                        j = j+1;
                    end
                    end_idx = j-1;
                    
                    %% Filter starts from here
                    % Condition 0-1: Discharge event length > 3 sec
                    I_seg = I(start_idx:end_idx);                    
                    if length(I_seg) < min_seg_length
                        i = end_idx+1;
                        continue;
                    end
                    
                    % Condition 0-2: Avg current > 1
                    if mean(abs(I_seg)) < idle_threshold
                        fprintf('[%s] Mean I %.2fA < 0.2A = 제외\n', dateStr, mean(abs(I_seg)));
                        i = end_idx+1;
                        continue;
                    end

                    % Condition 1: monotonous decrease, dIdt < 0, 단조감소에서 
                    % 전류가 구간내 양수가 있으면 제외로 변경
                    % if any(diff(I_seg) > 0) %  mean(diff(I_seg <= 0)) < 0.8  % 
                    if any(I_seg > 0)
                        fprintf('[%s] I > 0 구간 포함 제외\n', dateStr)
                        i = end_idx+1;
                        continue;
                    end

                    % Condition 2: Max dI: 0.1
                    if abs(max(I_seg) - min(I_seg)) > Cnom * 0.2 %   
                    % if any(abs(diff(I_seg))) > Cnom * 0.1
                        fprintf('[%s] ΔI %.2f > 0.1C = 제외\n', dateStr, abs(max(I_seg) - min(I_seg)));
                        i = end_idx+1;
                        continue;
                    end

                    % Condition 3: Max discharge current: 0.1C
                    if max(abs(I_seg)) > 1024*0.1
                        fprintf('[%s] Max I %.2fA > 0.1C = 제외\n', dateStr, max(abs(I_seg)));                       
                        i = end_idx+1;
                        continue
                    end
                    
                    % Condition 4: dSOC > 1%
                    % dSOC = abs(soc(end_idx) - soc(start_idx)); % max(soc(start_idx:end_idx) - min(soc(start_idx:end_idx)) < 1
                    dSOC = abs(max(soc) - min(soc));
                    if dSOC < min_dSOC % (1%)                    
                        fprintf('[%s] dSOC %.2f%% < %.2f%% = 제외\n', dateStr, dSOC, min_dSOC);
                        i = end_idx+1;
                        continue;
                    end

                    % Save 
                    eventCount = eventCount+1;
                    evtName = ['DchEvent_' , num2str(eventCount)];

                    Vcorrected.(dayField).(evtName).start_idx   = start_idx;
                    Vcorrected.(dayField).(evtName).end_idx     = end_idx;
                    Vcorrected.(dayField).(evtName).start_time  = t(start_idx);
                    Vcorrected.(dayField).(evtName).end_time    = t(end_idx);
                    Vcorrected.(dayField).(evtName).I_seq       = I(start_idx:end_idx);
                    Vcorrected.(dayField).(evtName).V_seq       = V(start_idx:end_idx);
                    Vcorrected.(dayField).(evtName).soc_seq     = soc(start_idx:end_idx);
                    Vcorrected.(dayField).(evtName).T_batt_seq  = T_batt(start_idx:end_idx);

                    fprintf('[%s] DchEvent %d Saved (%d seconds)\n', dateStr, eventCount, end_idx - start_idx + 1);
                    i = end_idx+1;
                else
                    i = i+1;
                end
            end
        end
    end
end

save(fullfile(saveDir, 'Vcorrected_struct.mat'), 'Vcorrected');

 
%% Discharge Event Plot
dateList = fieldnames(Vcorrected);
for i = 1:length(dateList)
    dayField = dateList{i};  % e.g., 'dischargeEvent_20230601'
    dateStr = extractAfter(dayField, 'dischargeEvent_');
    year = dateStr(1:4);
    month = dateStr(5:6);
    rawPath = fullfile(dataDir, year, [year month], ['Raw_' dateStr '.mat']);

    if ~exist(rawPath, 'file')
        fprintf('[%s] Raw 파일 없음 → 생략\n', dateStr);
        continue;
    end

    load(rawPath);
    t = Raw.sync_Time;
    I = Raw.Online_DC_Current;

    fig = figure('Name', ['Discharge Events - ', dateStr], 'Position', [100 100 1200 400]);
    plot(t, I, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
    ylabel('Current [A]'); xlabel('Time');
    title(['Raw Current with Discharge Events – ', dateStr]);
    grid on;

    events = fieldnames(Vcorrected.(dayField));
    for e = 1:length(events)
        event = Vcorrected.(dayField).(events{e});
        idx1 = event.start_idx;
        idx2 = event.end_idx;
        plot(t(idx1:idx2), I(idx1:idx2), 'r', 'LineWidth', 1.5);
    end

    % saveas(fig, fullfile(saveFigPath, ['DischargeEvents_' dateStr '.fig']));
    % close(fig);
end



%             fprintf('Processing %d discharge segments for %s\n', size(discharge_segments, 1), dateStr);
% 
%             % Handle if still in discharge at the end of data
%              if in_discharge && (length(I) - seg_start >= min_segment_length)
%                 discharge_segments(end+1,:) = [seg_start, length(I)];
%             end
% 
%             if isempty(discharge_segments)
%                 dischargeEventStruct.(safeDateStr) = struct([]);
%                 fprintf('[%s] No Discharge Event — Disregard Plot\n', dateStr);
%                 continue;
%             end
% 
%             % 방전 세그먼트와 정전류 방전 구간 시각화 추가
%             if ~isempty(discharge_segments)
%                 fig_discharge = figure('Name', ['Discharge Segments: ', dateStr], 'Position', [100, 100, 1600, 800]);
% 
%                 % 전류 subplot
%                 subplot(2,1,1);
%                 hold on; grid on;
%                 plot(t, I, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
%                 title(['Discharge Segments: ', dateStr]);
% 
%                 % 각 방전 세그먼트 하이라이트
%                 for s = 1:size(discharge_segments, 1)
%                     seg_start = discharge_segments(s, 1);
%                     seg_end = discharge_segments(s, 2);
%                     plot(t(seg_start:seg_end), I(seg_start:seg_end), 'r', 'LineWidth', 2);
%                 end
%                 ylabel('Current (A)');
% 
%                 % 전압 subplot
%                 subplot(2,1,2);
%                 hold on; grid on;
%                 plot(t, V, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
% 
%                 % 각 방전 세그먼트 하이라이트
%                 for s = 1:size(discharge_segments, 1)
%                     seg_start = discharge_segments(s, 1);
%                     seg_end = discharge_segments(s, 2);
%                     plot(t(seg_start:seg_end), V(seg_start:seg_end), 'r', 'LineWidth', 2);
%                 end
% 
%                 xlabel('Time');
%                 ylabel('Voltage (V)');
% 
%                 % 저장
%                 safeFigName = matlab.lang.makeValidName(['DischargeSegments_' dateStr]);
%                 saveas(fig_discharge, fullfile(saveDir, [safeFigName '.fig']));
%                 close(fig_discharge);
% 
%                 % 각 세그먼트의 상세 정보 출력
%                 fprintf('Discharge Segments Details:\n');
%                 for s = 1:size(discharge_segments, 1)
%                     seg_start = discharge_segments(s, 1);
%                     seg_end = discharge_segments(s, 2);
%                     seg_duration = seconds(t(seg_end) - t(seg_start));
%                     seg_mean_current = mean(I(seg_start:seg_end));
%                     seg_mean_voltage = mean(V(seg_start:seg_end));
% 
%                     fprintf('  Segment %d:\n', s);
%                     fprintf('    Start Time: %s\n', t(seg_start));
%                     fprintf('    End Time: %s\n', t(seg_end));
%                     fprintf('    Duration: %.2f seconds\n', seg_duration);
%                     fprintf('    Mean Current: %.2f A\n', seg_mean_current);
%                     fprintf('    Mean Voltage: %.2f V\n', seg_mean_voltage);
%                 end
%             end
% 
%             % 2. Process each discharge segment
%             for s = 1:size(discharge_segments, 1)
%                 start_idx = discharge_segments(s, 1);
%                 end_idx   = discharge_segments(s, 2);
% 
%                 % Extract segment data
%                 seg_I   = I(start_idx:end_idx);
%                 seg_V   = V(start_idx:end_idx);
%                 seg_SCO = soc(start_idx:end_idx);
%                 seg_T   = T_batt(start_idx:end_idx);
%                 seg_t   = t(start_idx:end_idx);
% 
%                 % Find current steps within segment
%                   step_indices = [];
%                 for j = 2:length(seg_I)
%                     if abs(seg_I(j) - seg_I(j-1)) >= min_deltaI
%                         step_indices = [step_indices, j];
%                     end
%                 end
%                 if s == 1
%                     max_di = 0;
%                     for j = 2:length(seg_I)
%                         max_di = max(max_di, abs(seg_I(j) - seg_I(j-1)));
%                     end
%                     fprintf('  First segment max dI: %.2f A (threshold: %.2f A)\n', max_di, min_deltaI);
%                 end
% 
%                 % Process each step for DCIR calculation
%                 for step_idx = 1:length(step_indices)
%                     st_idx = step_indices(step_idx);
% 
%                     % Need at leat min_segment length after step
%                     if st_idx + min_segment_length > length(seg_I)
%                         continue;
%                     end
% 
%                     % Check for constant current step after step
%                     cc_window = seg_I(st_idx:st_idx + min_segment_length - 1);
%                     cc_mean   = mean(cc_window);
%                     cc_var    = max(abs(cc_window - cc_mean) / abs(cc_mean));
% 
%                     if cc_var > cc_variation_limit
%                         fprintf('  Constant Current Window Mean: %.4f\n', cc_mean);
%                         fprintf('  CC Filter: Current variation %.4f > threshold %.4f\n', cc_var, cc_variation_limit);
%                         cc_filter_count = cc_filter_count + 1;
%                         continue;
%                     end
% 
%                     % Apply discharge event filter
%                     is_low = true;
%                     for j = st_idx:st_idx + min_segment_length - 1
%                         if j > 1 && abs(seg_I(j) - seg_I(j-1)) > dynamic_filter_threshold
%                             is_low = false;
%                             break;
%                         end
%                     end
% 
%                     if ~is_low
%                         fprintf('  Dynamic Filter: Large current change detected\n');
%                         dynamic_filter_count = dynamic_filter_count + 1;
%                         continue;
%                     end
% 
%                     thruput_window = st_idx: st_idx + min_segment_length - 1;
%                     time_diff = seconds(seg_t(thruput_window(end)) - seg_t(thruput_window(1)));
%                     thruput = abs(mean(seg_I(thruput_window)) * time_diff / 3600); % [Ah]
%                     relative_thruput = thruput / Cnom;
% 
%                     if relative_thruput < min_thruput
%                         fprintf('  Thruput Filter: Relative thruput %.4f < threshold %.4f\n', relative_thruput, min_thruput);
%                         thruput_filter_count = thruput_filter_count + 1;
%                         continue;
%                     end
% 
%                     % Calculate DCIR
%                     I1 = seg_I(st_idx - 1);
%                     I2 = seg_I(st_idx);
%                     V1 = seg_V(st_idx - 1);
%                     V2 = seg_V(st_idx);
% 
%                     dI = I2 - I1;
%                     dV = V2 - V1;
% 
%                     if abs(dI) < 0.00001
%                         deltaI_filter_count = deltaI_filter_count + 1;
%                         continue;
%                     end
% 
%                     valid_events = valid_events + 1;
% 
%                     DCIR = dV / dI;
%                     DCIR_mOhm = DCIR * 1000; % mOhm
% 
%                     current_SOC = seg_SOC(st_idx);
%                     current_t   = seg_T(st_idx);
% 
%                     soc_bin = max(1, min(floor(current_SOC) + 1 , length(SOC_bins) - 1));
%                     temp_bins = 10:5:45;
%                     [~, temp_bins] = min(abs(temp_bins - current_t));
% 
%                     yName = matlab.lang.makeValidName(year);
%                     if DCIR_mOhm > 0 && DCIR < 300
%                         DCR_lookup.(yName)(soc_bin, temp_bin) = DCIR_mOhm;
%                     end
% 
%                     V_corrected = seg_V - seg_I * DCIR;
% 
%                     % Store discharge events
%                     event_count = event_count + 1;
%                     dischargeEvents(event_count).start_idx = start_idx;
%                     dischargeEvents(event_count).end_idx = end_idx;
%                     % dischargeEvents(event_count).start_time = t(start_idx);
%                     % dischargeEvents(event_count).end_time = t(end_idx);
%                     dischargeEvents(event_count).DCIR_mOhm = DCIR_mOhm;
%                     dischargeEvents(event_count).Vocv = Vcorrected;
%                     dischargeEvents(event_count).SOC_seq = soc(start_idx:end_idx);
%                     dischargeEvents(event_count).I_seq = I(start_idx:end_idx);
%                     dischargeEvents(event_count).V_seq = V(start_idx:end_idx);
% 
%                     % Add to qOCV bins
%                     for k = 1:length(V_corrected)
%                         soc_val = seg_SOC(st_idx + k -1);
%                         v_val = V_corrected(k);
%                         bin_idx = find(soc_val >= SOC_bins(1:end-1) & soc_val < SOC_bins(2:end), 1);
%                         if ~isempty(bin_idx)
%                             bName = ['bin' num2str(SOC_bins(bin_idx))];
%                             qOCV_binned.(yName).(bName)(end+1) = v_val;                            
%                         end                    
%                     end
%                 end
%             end
% 
%             % Store the event data
%             if event_count > 0
%                 dischargeEventStruct.(safeDateStr) = dischargeEvents;
%                 fprintf('Processed %s - Found %d Discharge Events\n',dateStr, event_count);
%             else
%                 dischargeEventStruct.(safeDateStr) =struct([]);
%                 fprintf('[%s] No Valid Discharge Events Found\n',dateStr);                
%             end
%             fprintf('[%s] Summary: %d segments, %d steps, filters(CC=%d, Dynamic=%d, Thruput=%d, DeltaI=%d), Valid=%d\n', ...
%             dateStr, total_segments, total_steps, cc_filter_count, dynamic_filter_count, ...
%             thruput_filter_count, deltaI_filter_count, valid_events);
% 
%         end
%     end
% end

%% Generate Yearly qOCV-SOC Curve subplot

% fig1 = figure('Position', [50, 100, 1200, 300]);
% fig1 = figure('Position', [50, 100, 1200, 300]); %, 'WindowStyle', 'docked');
% years = fieldnames(qOCV_binned);
% if isempty(years)
%     warning('No valid qOCV data found for any year');
%     return;
% end
% SOC_mid = SOC_bins(1:end-1) + 0.5;
% colors = lines(length(years));
% 
% for k = 1:length(years)
%     y = years{k};
%     subplot(1, length(years), k); hold on; grid on;
% 
%     v_mean = zeros(size(SOC_mid));
%     v_std  = zeros(size(SOC_mid));
% 
%     for b = 1:length(SOC_mid)
%         bName = ['bin' num2str(SOC_bins(b))];
%         vals = qOCV_binned.(y).(bName);
%         if ~isempty(vals)
%             v_mean(b) = mean(vals);
%             v_std(b)  = std(vals);
%         else
%             v_mean(b) = NaN;
%             v_std(b)  = NaN;
%         end
%     end
% 
%     plot(SOC_mid, v_mean, '-o', 'Color', colors(k,:), 'LineWidth', 2);
%     patch([SOC_mid, fliplr(SOC_mid)], [v_mean+v_std, fliplr(v_mean-v_std)], colors(k,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
%     xlabel('SOC [%]');
%     xlim([0 100]);
%     xticks(0:10:100);
%     ylabel('V_{OCV} [V]');
%     ylim([840 950]);
%     yticks(840:10:960);
%     title(['qOCV - ', y]);
% end
% 
% sgtitle('Yearly qOCV - V_{OCV})');
% saveas(fig1, fullfile(saveDir, 'SOC_OCV_qOCV_subplot_by_year.fig'));
% close(fig1);
% 
% %% Generate DCR Heat Map
% 
% SOC_bins_plot = SOC_bins(1:end-1);
% Temp_bins = 15:5:35;     
% 
% for y = 1:length(yearList)
%     year = yearList{y};
%     yName = matlab.lang.makeValidName(year);
% 
%     if ~isfield(DCR_lookup, yName)
%         warning("No valid data for %s DCIR Heatmap", year);
%         continue;
%     end
% 
%     DCR_map = DCR_lookup.(year);
%     if all(isnan(DCR_map(:)))
%         warning("Not enough data for %s DCIR Heatmap", year);
%         continue;
%     end
% 
%     fig2 = figure('Position', [100 100 800 500]);
%     imagesc(SOC_bins_DCR(1:end-1), Temp_bins(1:end-1), DCR_map');
%     set(gca, 'YDir', 'normal');
%     xlabel('SOC [%]');
%     ylabel('Battery Temperature [°C]');
%     title(['DCIR Heatmap (SOC–Temp) – ', year]);
%     colormap(jet);
%     c = colorbar;
%     c.Label.String = 'DCIR [mΩ]';
%     colorbar;
% 
%     % DCR 범위 자동 설정
%     dcr_values = DCR_map(~isnan(DCR_map));
%     if ~isempty(dcr_values)
%         caxis([min(dcr_values) max(dcr_values)]);
%     end
% 
%     % 저장 경로
%     % fig_file = fullfile(saveDir, ['DCR_Heatmap_', year, '.fig']);
%     saveas(fig2, fullfile(saveDir, ['DCR_Heatmap_', year, '.fig']));    
%     % saveas(gcf, fig_file);
%     close(fig2);
% end
% 
% 
% %% Generate Events V, I Plot
% dateList = fieldnames(dischargeEventStruct);
% 
% for i = 1:length(dateList)
%     dateStr = dateList{i};
%     year = dateStr(2:5);  % 'x20230630' → '2023'
%     rawPath = fullfile(dataDir, year, dateStr(2:7), ['Raw_' dateStr(2:end) '.mat']);
% 
%     % if ~exist(rawPath, 'file')
%     %     fprintf('파일 없음: %s\n', rawPath);
%     %     continue;
%     % end
% 
%     load(rawPath);
%     t = Raw.sync_Time;
%     I = Raw.Online_DC_Current;
%     V = Raw.Total_Average_CV_Sum;
% 
%     events = dischargeEventStruct.(dateStr);
%     if ~isstruct(events) || isempty(events) % || ~isfield(events, 'start_idx')
%         fprintf('[%s] No Discharge Event — Disregard Plot \n', dateStr);
%         continue;
%     end
% 
%     % fig2 = figure('Name', ['Discharge Events ', dateStr], 'Position', [100, 100, 1200, 800]);
%     fig3 = figure('Name', ['Discharge Events ', dateStr], 'Position', [100, 100, 1200, 800]); %, 'WindowStyle', 'docked');
% 
%     % Current subplot
%     subplot(2,1,1); hold on; grid on;
%     % plot(t, I, 'k'); ylabel('Current [A]');
%     plot(t, I, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
%     title(['Current with Discharge Events: ', dateStr]);
% 
%     for e = 1:length(events)
%         sIdx = events(e).start_idx;
%         eIdx = events(e).end_idx;
%         plot(t(sIdx:eIdx), I(sIdx:eIdx), 'r', 'LineWidth', 1.5);
%     end
%     ylabel('Current (A)');
% 
%     % Voltage subplot
%     subplot(2,1,2); hold on; grid on;
%     % plot(t, V, 'k'); ylabel('Voltage [V]');
%     plot(t, V, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);  
%     title('Voltage with Discharge Events');
% 
%     for e = 1:length(events)
%         sIdx = events(e).start_idx;
%         eIdx = events(e).end_idx;
%         plot(t(sIdx:eIdx), V(sIdx:eIdx), 'r', 'LineWidth', 1.5);
% 
%         if isfield(events(e), 'Vcorr_seq')
%             plot(t(sIdx:eIdx), events(e).Vcorr_seq, 'b', 'LineWidth',1.5);
%         end
%     end
% 
%     xlabel('Time');
%     ylabel('Voltage (V)');
%     legend('Raw Voltage', 'Original', 'Corrected');
%     safeFigName = matlab.lang.makeValidName(['DischargeEvents_' dateStr]);
%     figPath = fullfile(saveDir, [safeFigName '.fig']);
% 
%     % disp(['저장 위치: ', figPath]);
%     saveas(fig3, fullfile(saveDir, ['DischargeEvents_', dateStr, '.fig']));
%     close(fig3);
%     fprintf('Saved %s (%d Discharge Events)\n', dateStr, length(events));
% end
% 
% fprintf('qOCV analysis completed and saved to %s\n', saveDir);
% 


%% Upload to Git
% current_dir = pwd;
% cd('ver01');
% cd('D:/JCW/KENTECH/Projects/KEPCO/ESS_Data_Preprocessing/qOCV_SOC');
% 
% system('git add ver01/');
% system('git commit -m "Add ver01 Folder');
% 
% % system('git remote add origin https://github.com/AutoVistion/qOCV_SOC.git');
% 
% system('git push -u origin main');
% 
% cd(current_dir);
% disp("Git Upload Complete");