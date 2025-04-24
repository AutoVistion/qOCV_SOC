%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOC- OCV Curve
% 30min, 1hr rest voltage = OCV
% Generate yearly SOC - OCV Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

%% Directory
clc; clear; close all;

% === 설정 ===
dataDir  = 'D:\JCW\KENTECH\Projects\KEPCO\ESS_Data_Preprocessing';
yearList = {'2021', '2022', '2023'};
saveDir  = fullfile(dataDir, 'ESS_Discharging_Events_DCR');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% dt_list = [3]; %[1, 3, 5, 10];
C_rate_limit       = 1024*0.1; % ICA, DVA에서는 0.05C가 일반적임
min_deltaI         = 5; 
min_segment_length = 10;  % C-rate limit 최소 구간 길이(초)
% dynamic_filter_threshold = 0.1 * 1024;  % C-rate의 10% (0.1C)



dischargeEventStruct = struct();
SOC_bins = 0:1:100;  % SOC binning

%% Yearly SOC bin
for y = 1:length(yearList)
    yName = matlab.lang.makeValidName(yearList{y});
    for b = 1:length(SOC_bins)-1
        bName = ['bin' num2str(SOC_bins(b))];
        qOCV_binned.(yName).(bName) = [];
    end
end

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

            load(matFilePath);
            t = Raw.sync_Time;
            I = Raw.Online_DC_Current;
            V = Raw.Total_Average_CV_Sum;
            soc = Raw.Total_Average_SOC;
            T_batt = Raw.Plc_Battery_Temperature_sync;

            % 이 부분부터 수정 시작 - 기존 코드 삭제
            dischargeEvents = [];
            event_count = 0;
            
            % 1. 먼저 전체 방전 구간 감지 (시각화용)
            i = 1;
            discharge_start = [];
            discharge_end = [];
            in_discharge = false;
            
            % 방전 구간 식별 (전류가 음수인 구간)
            for j = 1:length(I)
                if I(j) < 0 && ~in_discharge
                    in_discharge = true;
                    discharge_start = [discharge_start, j];
                elseif I(j) >= 0 && in_discharge
                    in_discharge = false;
                    discharge_end = [discharge_end, j-1];
                end
            end
            
            % 마지막 구간이 방전 중이면 끝점 추가
            if in_discharge
                discharge_end = [discharge_end, length(I)];
            end
            
            % 방전 구간이 없으면 다음 파일로
            if isempty(discharge_start)
                dischargeEventStruct.(safeDateStr) = struct([]);
                fprintf('[%s] No Discharge Event — Disregard Plot\n', dateStr);
                continue;
            end
            
            % 2. 각 방전 구간 내에서 저동적 세그먼트 식별 및 DCR 계산
            dynamic_filter_threshold = C_rate_limit * 2; % C-rate의 10% (0.1C)
            min_segment_length = 10; % 최소 세그먼트 길이
            
            for seg = 1:length(discharge_start)
                seg_start = discharge_start(seg);
                seg_end = discharge_end(seg);
                
                % 방전 구간 내 저동적 세그먼트 찾기
                i = seg_start;
                while i < seg_end
                    % C-rate 제한 확인
                    if abs(I(i)) >= C_rate_limit
                        i = i + 1;
                        continue;
                    end
                    
                    start_idx = i;
                    
                    % 저동적 세그먼트 끝점 찾기
                    j = i + 1;
                    while j <= seg_end && I(j) < 0 && abs(I(j)) < C_rate_limit
                        % 다이나믹 필터 적용
                        if abs(I(j) - I(j-1)) > dynamic_filter_threshold
                            break;
                        end
                        j = j + 1;
                    end
                    
                    end_idx = j - 1;
                    
                    % 최소 길이 확인
                    if (end_idx - start_idx + 1) >= min_segment_length
                        % DCR 계산
                        I1 = I(start_idx);
                        I2 = I(end_idx);
                        V1 = V(start_idx);
                        V2 = V(end_idx);
                        
                        deltaI = I2 - I1;
                        deltaV = V2 - V1;
                        
                        % 최소 전류 변화 확인
                        if abs(deltaI) >= min_deltaI
                            DCR = deltaV / deltaI;
                            DCR_mOhm = DCR * 1000;
                            
                            % Vcorr_seq 계산
                            Vcorr_seq = V(start_idx:end_idx) - I(start_idx:end_idx) * DCR;
                            
                            % 이벤트 저장
                            event_count = event_count + 1;
                            dischargeEvents(event_count).start_idx = start_idx;
                            dischargeEvents(event_count).end_idx = end_idx;
                            dischargeEvents(event_count).start_time = t(start_idx);
                            dischargeEvents(event_count).end_time = t(end_idx);
                            dischargeEvents(event_count).DCR_mOhm = DCR_mOhm;
                            dischargeEvents(event_count).Vcorr_seq = Vcorr_seq;
                            dischargeEvents(event_count).SOC_seq = soc(start_idx:end_idx);
                            dischargeEvents(event_count).I_seq = I(start_idx:end_idx);
                            dischargeEvents(event_count).V_seq = V(start_idx:end_idx);
                            








                            % qOCV bin 누적
                            SOC_seq = soc(start_idx:end_idx);
                            for k = 1:length(SOC_seq)
                                soc_val = SOC_seq(k);
                                v_val = Vcorr_seq(k);
                                for b = 1:length(SOC_bins)-1
                                    if soc_val >= SOC_bins(b) && soc_val < SOC_bins(b+1)
                                        bName = ['bin' num2str(SOC_bins(b))];
                                        yName = matlab.lang.makeValidName(year);
                                        qOCV_binned.(yName).(bName)(end+1) = v_val;
                                        break;
                                    end
                                end
                            end
                        end
                    end
                    
                    i = end_idx + 1;
                end
            end
            
            % 전체 방전 구간은 시각화를 위해 별도 저장 (옵션)
            for seg = 1:length(discharge_start)
                visualEvents(seg).start_idx = discharge_start(seg);
                visualEvents(seg).end_idx = discharge_end(seg);
            end

            if event_count > 0
                dischargeEventStruct.(safeDateStr) = dischargeEvents;
                fprintf('Processed %s - Found %d discharge events\n', dateStr, event_count);
            else
                dischargeEventStruct.(safeDateStr) = struct([]);
                fprintf('[%s] No Low-Dynamic Discharge Event — Disregard Plot\n', dateStr);
            end
            % 수정 끝
        end
    end
end

%% === 연도별 qOCV subplot 생성 ===
% fig1 = figure('Position', [50, 100, 1200, 300]);
fig1 = figure('Position', [50, 100, 1200, 300]); %, 'WindowStyle', 'docked');

years = fieldnames(qOCV_binned);
SOC_mid = SOC_bins(1:end-1) + 0.5;
colors = lines(length(years));

for k = 1:length(years)
    y = years{k};
    subplot(1, length(years), k); hold on; grid on;

    v_mean = zeros(size(SOC_mid));
    for b = 1:length(SOC_mid)
        bName = ['bin' num2str(SOC_bins(b))];
        vals = qOCV_binned.(y).(bName);
        if ~isempty(vals)
            v_mean(b) = mean(vals);
        else
            v_mean(b) = NaN;
        end
    end

    plot(SOC_mid, v_mean, '-o', 'Color', colors(k,:), 'LineWidth', 2);
    xlabel('SOC [%]');
    xlim([0 100]);
    xticks(0:10:100);
    ylabel('V_{OCV} [V]');
    ylim([840 950]);
    yticks(840:10:960);
    title(['qOCV - ', y]);
    xlim([0 100]);
end

sgtitle('연도별 qOCV (V_{corrected})');
saveas(fig1, fullfile(saveDir, 'SOC_OCV_qOCV_subplot_by_year.fig'));
close(fig1);

%% === 일자별 방전 이벤트 원시데이터 플롯 생성 ===
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
    if ~isstruct(events) || isempty(events) || ~isfield(events, 'start_idx')
        fprintf('[%s] No Discharge Event — Disregard Plot \n', dateStr);
        continue;
    end

    % fig2 = figure('Name', ['Discharge Events ', dateStr], 'Position', [100, 100, 1200, 800]);
    fig2 = figure('Name', ['Discharge Events ', dateStr], 'Position', [100, 100, 1200, 800]); %, 'WindowStyle', 'docked');

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

    % Voltage subplot
    subplot(2,1,2); hold on; grid on;
    % plot(t, V, 'k'); ylabel('Voltage [V]');
    plot(t, V, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);  % 원시 전압 회색
    title('Voltage with Discharge Events');

    for e = 1:length(events)
        sIdx = events(e).start_idx;
        eIdx = events(e).end_idx;
        plot(t(sIdx:eIdx), V(sIdx:eIdx), 'r', 'LineWidth', 1.5);
    end
    xlabel('Time');

    safeFigName = matlab.lang.makeValidName(['DischargeEvents_' dateStr]);
    figPath = fullfile(saveDir, [safeFigName '.fig']);

    % disp(['저장 위치: ', figPath]);
    saveas(fig2, figPath);
    close(fig2);
    fprintf('Saved %s (%d Discharge Events)\n', dateStr, length(events));
end

%% DCR 룩업 테이블 생성 및 시각화
% SOC 및 온도 구간 정의
SOC_bins_DCR = 0:5:100;  % 5% 단위
Temp_bins = 15:5:35;     % 5°C 단위

for y = 1:length(yearList)
    year = yearList{y};

    DCR_map = DCR_lookup.(year);
    if all(isnan(DCR_map(:)))
        warning("Not enough data for %s DCR heatmap", year);
        continue;
    end

    fig3 = figure('Position', [100 100 800 500]);
    imagesc(SOC_bins_DCR(1:end-1), Temp_bins(1:end-1), DCR_map');
    set(gca, 'YDir', 'normal');
    xlabel('SOC [%]');
    ylabel('Battery Temperature [°C]');
    title(['DCR Heatmap (SOC–Temp) – ', year]);
    colormap(jet);
    colorbar;

    % DCR 범위 자동 설정
    dcr_min = min(DCR_map(:), [], 'omitnan');
    dcr_max = max(DCR_map(:), [], 'omitnan');
    if ~isnan(dcr_min) && ~isnan(dcr_max)
        caxis([dcr_min dcr_max]);
    end

    % 저장 경로
    fig_file = fullfile(saveDir, ['DCR_Heatmap_', year, '.fig']);
    saveas(gcf, fig_file);
    close(gcf);
end
