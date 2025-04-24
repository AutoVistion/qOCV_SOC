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
yearList = {'2021'}; % , '2022', '2023'};
saveDir  = fullfile(dataDir, 'qOCV_SOC\ver01');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameters
Cnom                     = 1024;      % [Ah]
C_rate_limit             = 0.2*Cnom;  % 0.1C Limit for Discharge Events
min_deltaI               = 10;        % 0.05*Cnom; % minimum C-rate change
cc_variation_limit       = 0.65;      % for constant current detection(2%)
min_segment_length       = 3;         % C-rate limit 최소 구간 길이(초)
dynamic_filter_threshold = 0.2*Cnom;  % C-rate의 10% (0.1C)
min_thruput              = 0.001;     % charge/discharge thruput filter for SOC variation
SOC_bins                 = 0:1:100;   % SOC binning

dischargeEventStruct = struct();
qOCV_binned = struct();
DCR_lookup = struct();

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

            total_segments = 0;
            total_steps = 0;
            cc_filter_count = 0;
            dynamic_filter_count = 0;
            thruput_filter_count = 0;
            deltaI_filter_count = 0;
            valid_events = 0;

            load(matFilePath);
            t = Raw.sync_Time;
            I = Raw.Online_DC_Current;
            V = Raw.Total_Average_CV_Sum;
            soc = Raw.Total_Average_SOC;
            T_batt = Raw.Plc_Battery_Temperature_sync;

            dischargeEvents = [];
            event_count = 0;
                        
            % 1. Identify discharge segments
            discharge_segments = [];
            in_discharge = false;
            seg_start = 0;
            
            % Find constant current discharing phase
            for j = 1:length(I)
                if I(j) < 0 && ~in_discharge
                    in_discharge = true;
                    seg_start = j;
                elseif (I(j) >= 0 || j == length(I)) && in_discharge
                    in_discharge = false;
                    if j - seg_start >=min_segment_length
                        discharge_segments(end+1,:) = [seg_start, j-1];
                    end
                end
            end
            
            fprintf('Processing %d discharge segments for %s\n', size(discharge_segments, 1), dateStr);

            % Handle if still in discharge at the end of data
             if in_discharge && (length(I) - seg_start >= min_segment_length)
                discharge_segments(end+1,:) = [seg_start, length(I)];
            end
                        
            if isempty(discharge_segments)
                dischargeEventStruct.(safeDateStr) = struct([]);
                fprintf('[%s] No Discharge Event — Disregard Plot\n', dateStr);
                continue;
            end
              
            % 2. Process each discharge segment
            for s = 1:size(discharge_segments, 1)
                start_idx = discharge_segments(s, 1);
                end_idx   = discharge_segments(s, 2);
               
                % Extract segment data
                seg_I   = I(start_idx:end_idx);
                sev_V   = V(start_idx:end_idx);
                seg_SCO = soc(start_idx:end_idx);
                seg_T   = T_batt(start_idx:end_idx);
                seg_t   = t(start_idx:end_idx);

                % Find current steps within segment
                  step_indices = [];
                for j = 2:length(seg_I)
                    if abs(seg_I(j) - seg_I(j-1)) >= min_deltaI
                        step_indices = [step_indices, j];
                    end
                end
                if s == 1
                    max_di = 0;
                    for j = 2:length(seg_I)
                        max_di = max(max_di, abs(seg_I(j) - seg_I(j-1)));
                    end
                    fprintf('  First segment max dI: %.2f A (threshold: %.2f A)\n', max_di, min_deltaI);
                end

                if s == 1
                    max_di = 0;
                    for j = 2:length(seg_I)
                        max_di = max(max_di, abs(seg_I(j) - seg_I(j-1)));
                    end
                    fprintf('  First segment max dI: %.2f A (threshold: %.2f A)\n', max_di, min_deltaI);
                end

                % Process each step for DCIR calculation
                for step_idx = 1:length(step_indices)
                    st_idx = step_indices(step_idx);

                    % Need at leat min_segment length after step
                    if st_idx + min_segment_length > length(seg_I)
                        continue;
                    end

                    % Check for constant current step after step
                    cc_window = seg_I(st_idx:st_idx + min_segment_length - 1);
                    cc_mean   = mean(cc_window);
                    cc_var    = max(abs(cc_window - cc_mean) / abs(cc_mean));

                    if cc_var > cc_variation_limit
                        cc_filter_count = cc_filter_count + 1;
                        continue;
                    end
                    
                    % Apply discharge event filter
                    is_low = true;
                    for j = st_idx:st_idx + min_segment_length - 1
                        if j > 1 && abs(seg_I(j) - seg_I(j-1)) > dynamic_filter_threshold
                            is_low = false;
                            break;
                        end
                    end

                    if ~is_low
                        dynamic_filter_count = dynamic_filter_count + 1;
                        continue;
                    end

                    thruput_window = st_idx: st_idx + min_segment_length - 1;
                    thruput = abs(trapz(seg_t(thruput_window), seg_I(thruput_window))) / 3600; % [Ah]
                    relative_thruput = thruput / Cnom;

                    if relative_thruput < min_thruput
                        continue;
                    end

                    % Calculate DCIR
                    I1 = seg_I(st_idx - 1);
                    I2 = seg_I(st_idx);
                    V1 = seg_V(st_idx - 1);
                    V2 = seg_V(st_idx);

                    dI = I2 - I1;
                    dV = V2 - V1;

                    if abs(dI) < 0.00001
                        deltaI_filter_count = deltaI_filter_count + 1;
                        continue;
                    end
                    
                    valid_events = valid_events + 1;

                    DCIR = dV / dI;
                    DCIR_mOhm = DCIR * 1000; % mOhm

                    current_SOC = seg_SOC(st_idx);
                    current_t   = seg_T(st_idx);

                    soc_bin = max(1, min(floor(current_SOC) + 1 , length(SOC_bins) - 1));
                    temp_bins = 10:5:45;
                    [~, temp_bins] = min(abs(temp_bins - current_t));

                    yName = matlab.lang.makeValidName(year);
                    if DCIR_mOhm > 0 && DCIR < 300
                        DCR_lookup.(yName)(soc_bin, temp_bin) = DCIR_mOhm;
                    end

                    V_corrected = seg_V - seg_I * DCIR;

                    % Store discharge events
                    event_count = event_count + 1;
                    dischargeEvents(event_count).start_idx = start_idx;
                    dischargeEvents(event_count).end_idx = end_idx;
                    % dischargeEvents(event_count).start_time = t(start_idx);
                    % dischargeEvents(event_count).end_time = t(end_idx);
                    dischargeEvents(event_count).DCIR_mOhm = DCIR_mOhm;
                    dischargeEvents(event_count).Vocv = Vcorrected;
                    dischargeEvents(event_count).SOC_seq = soc(start_idx:end_idx);
                    dischargeEvents(event_count).I_seq = I(start_idx:end_idx);
                    dischargeEvents(event_count).V_seq = V(start_idx:end_idx);

                    % Add to qOCV bins
                    for k = 1:length(V_corrected)
                        soc_val = seg_SOC(st_idx + k -1);
                        v_val = V_corrected(k);
                        bin_idx = find(soc_val >= SOC_bins(1:end-1) & soc_val < SOC_bins(2:end), 1);
                        if ~isempty(bin_idx)
                            bName = ['bin' num2str(SOC_bins(bin_idx))];
                            qOCV_binned.(yName).(bName)(end+1) = v_val;                            
                        end                    
                    end
                end
            end

            % Store the event data
            if event_count > 0
                dischargeEventStruct.(safeDateStr) = dischargeEvents;
                fprintf('Processed %s - Found %d Discharge Events\n',dateStr, event_count);
            else
                dischargeEventStruct.(safeDateStr) =struct([]);
                fprintf('[%s] No Valid Discharge Events Found\n',dateStr);                
            end
            fprintf('[%s] Summary: %d segments, %d steps, filters(CC=%d, Dynamic=%d, Thruput=%d, DeltaI=%d), Valid=%d\n', ...
            dateStr, total_segments, total_steps, cc_filter_count, dynamic_filter_count, ...
            thruput_filter_count, deltaI_filter_count, valid_events);

        end
    end
end

%                  i = seg_start;
%                 while i < seg_end
%                     % C-rate 제한 확인
%                     if abs(I(i)) >= C_rate_limit
%                         i = i + 1;
%                         continue;
%                     end
% 
%                     start_idx = i;
% 
%                     % 저동적 세그먼트 끝점 찾기
%                     j = i + 1;
%                     while j <= seg_end && I(j) < 0 && abs(I(j)) < C_rate_limit
%                         % 다이나믹 필터 적용
%                         if abs(I(j) - I(j-1)) > dynamic_filter_threshold
%                             break;
%                         end
%                         j = j + 1;
%                     end
% 
%                     end_idx = j - 1;
% 
%                     % 최소 길이 확인
%                     if (end_idx - start_idx + 1) >= min_segment_length
%                         % DCR 계산
%                         I1 = I(start_idx);
%                         I2 = I(end_idx);
%                         V1 = V(start_idx);
%                         V2 = V(end_idx);
% 
%                         deltaI = I2 - I1;
%                         deltaV = V2 - V1;
% 
%                         % 최소 전류 변화 확인
%                         if abs(deltaI) >= min_deltaI
%                             DCR = deltaV / deltaI;
%                             DCR_mOhm = DCR * 1000;
% 
%                             % Vcorr_seq 계산
%                             Vcorr_seq = V(start_idx:end_idx) - I(start_idx:end_idx) * DCR;
% 
%                             % 이벤트 저장
%                             event_count = event_count + 1;
%                             dischargeEvents(event_count).start_idx = start_idx;
%                             dischargeEvents(event_count).end_idx = end_idx;
%                             dischargeEvents(event_count).start_time = t(start_idx);
%                             dischargeEvents(event_count).end_time = t(end_idx);
%                             dischargeEvents(event_count).DCR_mOhm = DCR_mOhm;
%                             dischargeEvents(event_count).Vcorr_seq = Vcorr_seq;
%                             dischargeEvents(event_count).SOC_seq = soc(start_idx:end_idx);
%                             dischargeEvents(event_count).I_seq = I(start_idx:end_idx);
%                             dischargeEvents(event_count).V_seq = V(start_idx:end_idx);
% 
%                             % qOCV bin 누적
%                             SOC_seq = soc(start_idx:end_idx);
%                             for k = 1:length(SOC_seq)
%                                 soc_val = SOC_seq(k);
%                                 v_val = Vcorr_seq(k);
%                                 for b = 1:length(SOC_bins)-1
%                                     if soc_val >= SOC_bins(b) && soc_val < SOC_bins(b+1)
%                                         bName = ['bin' num2str(SOC_bins(b))];
%                                         yName = matlab.lang.makeValidName(year);
%                                         qOCV_binned.(yName).(bName)(end+1) = v_val;
%                                         break;
%                                     end
%                                 end
%                             end
%                         end
%                     end
% 
%                     i = end_idx + 1;
%                 end
%             end
% 
%             % 전체 방전 구간은 시각화를 위해 별도 저장 (옵션)
%             for seg = 1:length(discharge_segments)
%                 visualEvents(seg).start_idx = discharge_segments(seg);
%                 visualEvents(seg).end_idx = discharge_end(seg);
%             end
% 
%             if event_count > 0
%                 dischargeEventStruct.(safeDateStr) = dischargeEvents;
%                 fprintf('Processed %s - Found %d discharge events\n', dateStr, event_count);
%             else
%                 dischargeEventStruct.(safeDateStr) = struct([]);
%                 fprintf('[%s] No Low-Dynamic Discharge Event — Disregard Plot\n', dateStr);
%             end
%             % 수정 끝
%         end
%     end
% end

%% Generate Yearly qOCV-SOC Curve subplot

% fig1 = figure('Position', [50, 100, 1200, 300]);
fig1 = figure('Position', [50, 100, 1200, 300]); %, 'WindowStyle', 'docked');
years = fieldnames(qOCV_binned);
if isempty(years)
    warning('No valid qOCV data found for any year');
    return;
end
SOC_mid = SOC_bins(1:end-1) + 0.5;
colors = lines(length(years));

for k = 1:length(years)
    y = years{k};
    subplot(1, length(years), k); hold on; grid on;

    v_mean = zeros(size(SOC_mid));
    v_std  = zeros(size(SOC_mid));

    for b = 1:length(SOC_mid)
        bName = ['bin' num2str(SOC_bins(b))];
        vals = qOCV_binned.(y).(bName);
        if ~isempty(vals)
            v_mean(b) = mean(vals);
            v_std(b)  = std(vals);
        else
            v_mean(b) = NaN;
            v_std(b)  = NaN;
        end
    end

    plot(SOC_mid, v_mean, '-o', 'Color', colors(k,:), 'LineWidth', 2);
    patch([SOC_mid, fliplr(SOC_mid)], [v_mean+v_std, fliplr(v_mean-v_std)], colors(k,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    xlabel('SOC [%]');
    xlim([0 100]);
    xticks(0:10:100);
    ylabel('V_{OCV} [V]');
    ylim([840 950]);
    yticks(840:10:960);
    title(['qOCV - ', y]);
end

sgtitle('Yearly qOCV - V_{OCV})');
saveas(fig1, fullfile(saveDir, 'SOC_OCV_qOCV_subplot_by_year.fig'));
close(fig1);

%% Generate DCR Heat Map

SOC_bins_plot = SOC_bins(1:end-1);
Temp_bins = 15:5:35;     

for y = 1:length(yearList)
    year = yearList{y};
    yName = matlab.lang.makeValidName(year);

    if ~isfield(DCR_lookup, yName)
        warning("No valid data for %s DCIR Heatmap", year);
        continue;
    end

    DCR_map = DCR_lookup.(year);
    if all(isnan(DCR_map(:)))
        warning("Not enough data for %s DCIR Heatmap", year);
        continue;
    end

    fig2 = figure('Position', [100 100 800 500]);
    imagesc(SOC_bins_DCR(1:end-1), Temp_bins(1:end-1), DCR_map');
    set(gca, 'YDir', 'normal');
    xlabel('SOC [%]');
    ylabel('Battery Temperature [°C]');
    title(['DCIR Heatmap (SOC–Temp) – ', year]);
    colormap(jet);
    c = colorbar;
    c.Label.String = 'DCIR [mΩ]';
    colorbar;

    % DCR 범위 자동 설정
    dcr_values = DCR_map(~isnan(DCR_map));
    if ~isempty(dcr_values)
        caxis([min(dcr_values) max(dcr_values)]);
    end

    % 저장 경로
    % fig_file = fullfile(saveDir, ['DCR_Heatmap_', year, '.fig']);
    saveas(fig2, fullfile(saveDir, ['DCR_Heatmap_', year, '.fig']));    
    % saveas(gcf, fig_file);
    close(fig2);
end


%% Generate Events V, I Plot
dateList = fieldnames(dischargeEventStruct);

for i = 1:length(dateList)
    dateStr = dateList{i};
    year = dateStr(2:5);  % 'x20230630' → '2023'
    rawPath = fullfile(dataDir, year, dateStr(2:7), ['Raw_' dateStr(2:end) '.mat']);

    % if ~exist(rawPath, 'file')
    %     fprintf('파일 없음: %s\n', rawPath);
    %     continue;
    % end

    load(rawPath);
    t = Raw.sync_Time;
    I = Raw.Online_DC_Current;
    V = Raw.Total_Average_CV_Sum;

    events = dischargeEventStruct.(dateStr);
    if ~isstruct(events) || isempty(events) % || ~isfield(events, 'start_idx')
        fprintf('[%s] No Discharge Event — Disregard Plot \n', dateStr);
        continue;
    end

    % fig2 = figure('Name', ['Discharge Events ', dateStr], 'Position', [100, 100, 1200, 800]);
    fig3 = figure('Name', ['Discharge Events ', dateStr], 'Position', [100, 100, 1200, 800]); %, 'WindowStyle', 'docked');

    % Current subplot
    subplot(2,1,1); hold on; grid on;
    % plot(t, I, 'k'); ylabel('Current [A]');
    plot(t, I, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    title(['Current with Discharge Events: ', dateStr]);

    for e = 1:length(events)
        sIdx = events(e).start_idx;
        eIdx = events(e).end_idx;
        plot(t(sIdx:eIdx), I(sIdx:eIdx), 'r', 'LineWidth', 1.5);
    end
    ylabel('Current (A)');

    % Voltage subplot
    subplot(2,1,2); hold on; grid on;
    % plot(t, V, 'k'); ylabel('Voltage [V]');
    plot(t, V, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);  
    title('Voltage with Discharge Events');

    for e = 1:length(events)
        sIdx = events(e).start_idx;
        eIdx = events(e).end_idx;
        plot(t(sIdx:eIdx), V(sIdx:eIdx), 'r', 'LineWidth', 1.5);

        if isfield(events(e), 'Vcorr_seq')
            plot(t(sIdx:eIdx), events(e).Vcorr_seq, 'b', 'LineWidth',1.5);
        end
    end

    xlabel('Time');
    ylabel('Voltage (V)');
    legend('Raw Voltage', 'Original', 'Corrected');
    safeFigName = matlab.lang.makeValidName(['DischargeEvents_' dateStr]);
    figPath = fullfile(saveDir, [safeFigName '.fig']);

    % disp(['저장 위치: ', figPath]);
    saveas(fig3, fullfile(saveDir, ['DischargeEvents_', dateStr, '.fig']));
    close(fig3);
    fprintf('Saved %s (%d Discharge Events)\n', dateStr, length(events));
end

fprintf('qOCV analysis completed and saved to %s\n', saveDir);



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