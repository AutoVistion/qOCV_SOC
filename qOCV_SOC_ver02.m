%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates SOC - qOCV Curve from Raw Data.
% 1. Find discharging event within 0.1C with constant current phase duration
% over 3 seconds.
% 2. Calculate DCIR from discharge events (using 3-second point)
% 3. Generate SOC qOCV Curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Directory Setup
dataDir  = 'D:\JCW\KENTECH\Projects\KEPCO\ESS_Data_Preprocessing';
yearList = {'2021'};
saveDir  = fullfile(dataDir, 'qOCV_SOC\ver01');

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameter Setting
Cnom            = 1024;      % [Ah]
idle_threshold  = 1;         % Current threshold for idle state [A]
min_seg_length  = 5;        % Minimum discharge length [seconds]
min_dSOC        = 1;         % Minimum SOC change [%]
% soc_bins        = 0:1:100;  % SOC binning


all_events = struct();

%% Main Processing Loop
fprintf('Processing Data...\n');
fprintf('=================\n');

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
            fieldName = ['date_' dateStr];

            % Load data
            % try
            load(matFilePath);
            t = Raw.sync_Time;
            I = Raw.Online_DC_Current;
            V = Raw.Total_Average_CV_Sum;
            soc = Raw.Total_Average_SOC;
            T_batt = Raw.Plc_Battery_Temperature_sync;

            % Detect discharge events
            discharge_phases = detectDischargePhases(t, I, V, soc, T_batt, min_seg_length, idle_threshold, Cnom, min_dSOC);

            % raw 데이터만
            all_events.(fieldName) = struct();
            phases = fieldnames(discharge_phases);
            for p = 1:length(phases)
                phase = discharge_phases.(phases{p});
                all_events.(fieldName).(phases{p}) = struct();
                all_events.(fieldName).(phases{p}).I_seq = phase.I_seq;
                all_events.(fieldName).(phases{p}).V_seq = phase.V_seq;
                all_events.(fieldName).(phases{p}).soc_seq = phase.soc_seq;
                all_events.(fieldName).(phases{p}).T_batt_seq = phase.T_batt_seq;
                all_events.(fieldName).(phases{p}).start_idx = phase.start_idx;
                all_events.(fieldName).(phases{p}).end_idx = phase.end_idx;
                all_events.(fieldName).(phases{p}).start_time = phase.start_time;
                all_events.(fieldName).(phases{p}).end_time = phase.end_time;
                all_events.(fieldName).(phases{p}).DCIR = phase.DCIR;
                all_events.(fieldName).(phases{p}).Vcorr_seq = phase.Vcorr_seq;
            end
        end
    end
end

save(fullfile(saveDir, 'all_events.mat'), 'all_events');

%% Generate Plots
fprintf('\nGenerating Plots...\n');
fprintf('=================\n');

disp(fieldnames(all_events))

dateList = fieldnames(all_events);
all_soc = [];
all_vcorr = [];
for d = 1:length(dateList)
    date = dateList{d};
    events = fieldnames(all_events.(date));
    for e = 1:length(events)
        evt = all_events.(date).(events{e});
        if isfield(evt, 'Vcorr_seq') && ~isempty(evt.Vcorr_seq)
            all_soc = [all_soc; evt.soc_seq(:)];
            all_vcorr = [all_vcorr; evt.Vcorr_seq(:)];
        end
    end
end

plotValidQOCVCurves(all_events, saveDir);

% plotDischargeEventsOverlay_ver02(dataDir, all_events, saveDir)

plotMonthlyQOCVCurves(all_events, saveDir)

plotDCIR_SOC_Tbatt_heatmap_sparse(all_events, saveDir)

% disp([soc_list tbatt_list dcir_list])
% figure; scatter(soc_list, tbatt_list, 50, dcir_list, 'filled');
% xlabel('SOC'); ylabel('Tbatt'); colorbar; title('이벤트별 SOC-Tbatt 분포');

% disp(length(soc_list))

% tabulate(soc_idx)
% tabulate(tbatt_idx)

fprintf('\nAnalysis Complete\n');

%% Function Definitions

function discharge_phases = detectDischargePhases(t, I, V, soc, T_batt, min_seg_length, idle_threshold, Cnom, min_dSOC)
    discharge_phases = struct();
    eventCount = 0;
    min_discharge_current = 2;
    current_jump = 30; % 급격한 변화 임계값 (A)
    i = 2;
    while i <= length(I) - min_seg_length
        if I(i) < 0 && abs(I(i-1)) <= idle_threshold
            fprintf('[%s] 방전이벤트 후보: i = %d, I = %.2f, I(i-1) = %.2f\n', datestr(t(i), 'yyyymmdd'), i, I(i), I(i-1));
            k = i - 1;
            while k >= 1 && abs(I(k)) <= idle_threshold
                k = k - 1;
            end
            idle_start = k + 1;
            idle_range = idle_start:i-1;

            if ~isempty(idle_range) && idle_range(end) < length(I)
                % idle 구간 마지막 값
                last_idle_idx = idle_range(end);
                % 바로 다음 값이 min_discharge_current 이상(음수 방향)인지 확인
                if I(last_idle_idx + 1) <= -min_discharge_current
                    start_idx = last_idle_idx;
                else
                    % 조건 불충족 시, 이벤트 후보로 인정하지 않음 (continue 등)
                    i = i + 1;
                    continue;
                end
            else
                start_idx = i; % fallback
            end

            % 종료조건: 전류가 0 이상이 되거나, dI > current_jump(양수 방향)
            j = i+1;
            while j <= length(I)
                dI = I(j) - I(j-1);
                if I(j) >= 0 || dI > current_jump
                    end_idx = j-1;
                    break;
                end
                
                j = j+1;
            end
            % end_idx = j-1;

            %% Filter starts from here
            % Condition 0-1: Discharge event length > 3 sec
            I_seg = I(start_idx:end_idx);
            if length(I_seg) < min_seg_length
                i = end_idx+1;
                continue;
            end

            % Condition 0-2: Avg current > min_discharge_current
            if mean(abs(I_seg)) < min_discharge_current
                fprintf('[%s] Mean I %.2fA < %.2fA = 제외\n', datestr(t(i), 'yyyymmdd'), mean(abs(I_seg)), min_discharge_current);
                i = end_idx+1;
                continue;
            end

            % Condition 0-3: Max current > min_discharge_current
            if max(abs(I_seg)) < min_discharge_current
                fprintf('[%s] Max I %.2fA < %.2fA = 제외\n', datestr(t(i), 'yyyymmdd'), max(abs(I_seg)), min_discharge_current);
                i = end_idx+1;
                continue;
            end

            % Condition 1: monotonous decrease, dIdt < 0, 단조감소에서
            % 전류가 구간내 양수가 있으면 제외로 변경
            if any(I_seg > 0)
                fprintf('[%s] I > 0 구간 포함 제외\n', datestr(t(i), 'yyyymmdd'))
                i = end_idx+1;
                continue;
            end

            % % Condition 2: Max dI: 0.1
            % if abs(max(I_seg) - min(I_seg)) > Cnom * 0.1
            %     fprintf('[%s] ΔI %.2f > 0.1C = 제외\n', datestr(t(i), 'yyyymmdd'), abs(max(I_seg) - min(I_seg)));
            %     i = end_idx+1;
            %     continue;
            % end

            % Condition 3: Max discharge current: 0.1C
            if max(abs(I_seg)) > Cnom * 0.25
                fprintf('[%s] Max I %.2fA > 0.1C = 제외\n', datestr(t(i), 'yyyymmdd'), max(abs(I_seg)));
                i = end_idx+1;
                continue
            end

            % Condition 4: dSOC > min_dSOC
            soc_seg = soc(start_idx:end_idx);
            dSOC = abs(soc_seg(end) - soc_seg(1));
            if dSOC < min_dSOC
                fprintf('[%s] dSOC %.2f%% < %.2f%% = 제외\n', datestr(t(i), 'yyyymmdd'), dSOC, min_dSOC);
                i = end_idx+1;
                continue;
            end

            % Save
            eventCount = eventCount+1;
            evtName = ['DchEvent_', num2str(eventCount)];

            % DCIR 및 Vcorr 계산 추가
            I_evt = I(start_idx:end_idx);
            V_evt = V(start_idx:end_idx);
            V1 = V_evt(1); I1 = I_evt(1);
            V2 = V_evt(end); I2 = I_evt(end);
            dI = I2 - I1;
            if abs(dI) > 1e-6
                DCIR = (V2 - V1) / dI;
                Vcorr_seq = V_evt - I_evt * DCIR;
            else
                DCIR = NaN;
                Vcorr_seq = nan(size(V_evt));
            end

            discharge_phases.(evtName).start_idx = start_idx;
            discharge_phases.(evtName).end_idx = end_idx;
            discharge_phases.(evtName).start_time = t(start_idx);
            discharge_phases.(evtName).end_time = t(end_idx);
            discharge_phases.(evtName).I_seq = I_evt;
            discharge_phases.(evtName).V_seq = V_evt;
            discharge_phases.(evtName).soc_seq = soc(start_idx:end_idx);
            discharge_phases.(evtName).T_batt_seq = T_batt(start_idx:end_idx);
            discharge_phases.(evtName).DCIR = DCIR;
            discharge_phases.(evtName).Vcorr_seq = Vcorr_seq;

            fprintf('[%s] DchEvent %d Saved (%d seconds)\n', datestr(t(i), 'yyyymmdd'), eventCount, end_idx - start_idx + 1);
            i = end_idx+1;
        else
            i = i+1;
        end
    end
end

%%
function plotValidQOCVCurves(all_events, saveDir)
    % 1. 각 이벤트별 부분 qOCV 곡선(shift 전)
    % 2. 임시 평균 qOCV 곡선(shift 전)
    % 3. shift 후 최종 평균 qOCV 곡선(논문 방식)

    all_soc = [];
    all_vcorr = [];
    event_soc = {};
    event_vcorr = {};
    dateList = fieldnames(all_events);
    colors = lines(100);
    event_count = 0;
    for d = 1:length(dateList)
        events = fieldnames(all_events.(dateList{d}));
        for e = 1:length(events)
            evt = all_events.(dateList{d}).(events{e});
            if isfield(evt, 'Vcorr_seq') && ~isempty(evt.Vcorr_seq)
                all_soc = [all_soc; evt.soc_seq(:)];
                all_vcorr = [all_vcorr; evt.Vcorr_seq(:)];
                event_count = event_count + 1;
                event_soc{event_count} = evt.soc_seq(:);
                event_vcorr{event_count} = evt.Vcorr_seq(:);
            end
        end
    end
    if isempty(all_soc) || isempty(all_vcorr)
        warning('No valid qOCV data to plot.');
        return;
    end

    %% 1. 각 이벤트별 부분 qOCV 곡선(shift 전)
    fig1 = figure('Position', [100, 100, 800, 600]); hold on; grid on;
    for k = 1:event_count
        soc_seq = event_soc{k};
        plot(soc_seq, event_vcorr{k}, '-', 'Color', colors(mod(k-1, size(colors,1))+1,:), 'LineWidth', 1.5);
    end
    xlabel('SOC [%]'); ylabel('qOCV [V]');
    title('Partial qOCV Curves (All Events, Before Alignment)');
    xlim([0 100]);
    saveas(fig1, fullfile(saveDir, 'Partial_qOCV_Curves_BeforeAlignment.fig'));
    close(fig1);

    %% 2. 임시 평균 qOCV 곡선(shift 전)
    all_soc = [];
    for k = 1:event_count
        all_soc = [all_soc; event_soc{k}(:)];
    end
    soc_bins = 0:1:100;
    mean_v = nan(size(soc_bins));
    for i = 1:length(soc_bins)
        idx = all_soc >= soc_bins(i) & all_soc < soc_bins(i)+1;
        if sum(idx) >= 1
            mean_v(i) = mean(all_vcorr(idx));
        end
    end
    fig2 = figure('Position', [100, 100, 800, 600]); hold on; grid on;
    scatter(all_soc, all_vcorr, 5, [0.7 0.7 0.7], 'filled');
    plot(soc_bins, mean_v, 'b-', 'LineWidth', 2);
    xlabel('SOC [%]'); ylabel('qOCV [V]');
    title('Initial Mean qOCV Curve (Before Alignment)');
    xlim([0 100]);
    legend('All Points','Mean qOCV','Location','best');
    saveas(fig2, fullfile(saveDir, 'Mean_qOCV_Curve_BeforeAlignment.fig'));
    close(fig2);

    aligned_soc = [];
    aligned_vcorr = [];
    for k = 1:event_count
        soc_k = event_soc{k};
        vcorr_k = event_vcorr{k};
        shifts = -5:0.1:5;
        best_shift = 0;
        min_err = inf;
        for s = 1:length(shifts)
            soc_shifted = soc_k + shifts(s);
            % 평균 곡선과 겹치는 구간만 비교
            idx_overlap = soc_shifted >= min(soc_bins) & soc_shifted <= max(soc_bins);
            soc_overlap = soc_shifted(idx_overlap);
            vcorr_overlap = vcorr_k(idx_overlap);
            % 평균 곡선에서 해당 SOC에 대응하는 값 추출
            mean_v_interp = interp1(soc_bins, mean_v, soc_overlap, 'linear', 'extrap');
            % 오차 계산 (MAE)
            err = mean(abs(vcorr_overlap - mean_v_interp), 'omitnan');
            if err < min_err
                min_err = err;
                best_shift = shifts(s);
            end
        end
        soc_aligned = soc_k + best_shift;
        aligned_soc = [aligned_soc; soc_aligned(:)];
        aligned_vcorr = [aligned_vcorr; vcorr_k(:)];
    end
    mean_v_aligned = nan(size(soc_bins));
    for i = 1:length(soc_bins)
        idx = aligned_soc >= soc_bins(i) & aligned_soc < soc_bins(i)+1;
        if sum(idx) >= 1
            mean_v_aligned(i) = mean(aligned_vcorr(idx));
        end
    end
    fig3 = figure('Position', [100, 100, 800, 600]); hold on; grid on;
    scatter(aligned_soc, aligned_vcorr, 5, [0.7 0.7 0.7], 'filled');
    plot(soc_bins, mean_v_aligned, 'r-', 'LineWidth', 2);
    xlabel('SOC [%]'); ylabel('qOCV [V]');
    title('Final Mean qOCV Curve (After Alignment, Paper Method)');
    xlim([0 100]);
    legend('All Points','Mean qOCV (Aligned)','Location','best');
    saveas(fig3, fullfile(saveDir, 'Mean_qOCV_Curve_AfterAlignment.fig'));
    close(fig3);
end

%%
% function plotDischargeEventsOverlay_ver02(dataDir, all_events, saveDir)
% 
%     dateList = fieldnames(all_events);
% 
%     for d = 1:length(dateList)
%         dayField = dateList{d}; % 
%         dateStr = extractAfter(dayField, 'date_'); % '20210101'
%         year = dateStr(1:4);
%         month = dateStr(5:6);
%         rawPath = fullfile(dataDir, year, [year month], ['Raw_' dateStr '.mat']);
% 
%         if ~exist(rawPath, 'file')
%             fprintf('[%s] Raw file not found\n', dateStr);
%             continue;
%         end
% 
%         load(rawPath); % Raw 구조체 포함
%         t = Raw.sync_Time;
%         I = Raw.Online_DC_Current;
%         V = Raw.Total_Average_CV_Sum;
%         soc = Raw.Total_Average_SOC;
% 
%         fig = figure('Visible', 'on', 'Name', ['Discharge Events Overlay - ', dateStr], 'Position', [100 100 1400 900]);
% 
%         % 1. Current subplot
%         subplot(3,1,1);
%         plot(t, I, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
%         ylabel('Current [A]'); xlabel('Time');
%         title(['Raw Current with Detected Discharge Events – ', dateStr]);
%         grid on;
% 
%         % 2. Voltage subplot
%         subplot(3,1,2);
%         plot(t, V, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
%         ylabel('Voltage [V]'); xlabel('Time');
%         title(['Raw Voltage with Detected Discharge Events – ', dateStr]);
%         grid on;
% 
%         % 3. SOC subplot
%         subplot(3,1,3);
%         plot(t, soc, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
%         ylabel('SOC [%]'); xlabel('Time');
%         title(['Raw SOC with Detected Discharge Events – ', dateStr]);
%         grid on;
% 
%         % Overlay detected events in red
%         eventNames = fieldnames(all_events.(dayField));
%         for e = 1:length(eventNames)
%             evt = all_events.(dayField).(eventNames{e});
%             idx1 = evt.start_idx;
%             idx2 = evt.end_idx;
% 
%             subplot(3,1,1);
%             plot(t(idx1:idx2), I(idx1:idx2), 'r', 'LineWidth', 2);
% 
%             subplot(3,1,2);
%             plot(t(idx1:idx2), V(idx1:idx2), 'r', 'LineWidth', 2);
% 
%             subplot(3,1,3);
%             plot(t(idx1:idx2), soc(idx1:idx2), 'r', 'LineWidth', 2);
%         end
% 
%         % Save figure
%         saveas(fig, fullfile(saveDir, ['DischargeEventsOverlay_' dateStr '.fig']));
%         close(fig);
%     end
% end

function plotMonthlyQOCVCurves(all_events, saveDir)
    % 월별로 공칭용량 추정 및 qOCV-SOCnom 곡선 플롯
    dateList = fieldnames(all_events);
    months = strings(1, length(dateList));
    for d = 1:length(dateList)
        date = extractAfter(dateList{d}, 'date_');
        months(d) = date(1:6); % YYYYMM
    end
    unique_months = unique(months);
    colors = lines(length(unique_months));
    fig = figure('Position', [100, 100, 800, 600]); hold on; grid on;
    legends = {};
    for m = 1:length(unique_months)
        % 해당 월의 이벤트 모으기
        idx_month = months == unique_months(m);
        socnom_all = [];
        vcorr_all = [];
        cnom_list = [];
        for d = find(idx_month)
            events = fieldnames(all_events.(dateList{d}));
            for e = 1:length(events)
                evt = all_events.(dateList{d}).(events{e});
                if isfield(evt, 'Vcorr_seq') && ~isempty(evt.Vcorr_seq)
                    I_seq = evt.I_seq(:);
                    Vcorr_seq = evt.Vcorr_seq(:);
                    % 누적 방전량 계산 (Ah)
                    dAh_seq = cumsum(I_seq)/3600; % 샘플링 간격 1초 가정
                    dAh_evt = abs(trapz(I_seq)/3600);
                    dSOC = abs(evt.soc_seq(end) - evt.soc_seq(1));
                    if dSOC > 0.1
                        cnom_evt = dAh_evt / (dSOC/100);
                        cnom_list = [cnom_list; cnom_evt];
                    end
                end
            end
        end
        if isempty(cnom_list)
            continue;
        end
        % 월별 평균 Cnom
        cnom_month = mean(cnom_list, 'omitnan');
        for d = find(idx_month)
            events = fieldnames(all_events.(dateList{d}));
            for e = 1:length(events)
                evt = all_events.(dateList{d}).(events{e});
                if isfield(evt, 'Vcorr_seq') && ~isempty(evt.Vcorr_seq)
                    I_seq = evt.I_seq(:);
                    Vcorr_seq = evt.Vcorr_seq(:);
                    dAh_seq = cumsum(I_seq)/3600; % 누적 방전량
                    SOCnom_seq = 100 - 100 * abs(dAh_seq) / cnom_month;

                    SOCnom_seq_shifted = SOCnom_seq - SOCnom_seq(1) + evt.soc_seq(1);
                    socnom_all = [socnom_all; SOCnom_seq_shifted];
                    vcorr_all = [vcorr_all; Vcorr_seq];
                end
            end
        end
        if isempty(socnom_all) || isempty(vcorr_all)
            continue;
        end
        % SOCnom binning, 평균 qOCV
        soc_bins = 0:1:100;
        mean_v = nan(size(soc_bins));
        for i = 1:length(soc_bins)
            idx = socnom_all >= soc_bins(i) & socnom_all < soc_bins(i)+1;
            if sum(idx) >= 1
                mean_v(i) = mean(vcorr_all(idx));
            end
        end
        plot(soc_bins, mean_v, '-o', 'Color', colors(m,:), 'LineWidth', 2, 'MarkerSize', 4);
        legends{end+1} = sprintf('%s (Cnom=%.0fAh)', unique_months(m), cnom_month);
    end
    xlabel('SOC_{Nominal} in %'); ylabel('qOCV in V');
    title('Monthly Mean qOCV Curves with Estimated C_{nom} (True SOCnom)');
    xlim([0 100]);
    ylim auto;
    legend(legends, 'Location', 'best');
    saveas(fig, fullfile(saveDir, 'Monthly_qOCV_Curves_with_Cnom.fig'));
    close(fig);
end

%%
function plotDCIR_SOC_Tbatt_heatmap_sparse(all_events, saveDir)
    % 2D 히트맵: X축 SOC, Y축 Tbatt, 값은 DCIR(mΩ) (실제 데이터 있는 구간만 표시)
    dateList = fieldnames(all_events);
    soc_list = [];
    tbatt_list = [];
    dcir_list = [];
    for d = 1:length(dateList)
        events = fieldnames(all_events.(dateList{d}));
        for e = 1:length(events)
            evt = all_events.(dateList{d}).(events{e});
            if isfield(evt, 'DCIR') && ~isnan(evt.DCIR) && isfield(evt, 'soc_seq') && isfield(evt, 'T_batt_seq')
                soc_val = evt.soc_seq(1); % 이벤트 시작 SOC
                tbatt_val = mean(evt.T_batt_seq); % 이벤트 평균 Tbatt
                dcir_val = evt.DCIR * 1000; % [Ohm] -> [mΩ]
                soc_list = [soc_list; soc_val];
                tbatt_list = [tbatt_list; tbatt_val];
                dcir_list = [dcir_list; dcir_val];
            end
        end
    end

    % 원하는 bin 설정
    soc_bins = 0:10:100;
    tbatt_bins = 15:1:40;

    % 각 이벤트의 SOC, Tbatt가 어느 bin에 속하는지 인덱스 계산
    [~, soc_idx] = histc(soc_list, soc_bins);
    [~, tbatt_idx] = histc(tbatt_list, tbatt_bins);

    % 실제 데이터가 있는 bin만 추출
    bin_table = table(soc_idx, tbatt_idx, dcir_list, 'VariableNames', {'SOC', 'Tbatt', 'DCIR'});
    bin_table = bin_table(soc_idx > 0 & tbatt_idx > 0, :);

    % 각 bin별 평균 DCIR 계산
    [unique_bins, ~, ic] = unique([bin_table.SOC, bin_table.Tbatt], 'rows');
    mean_dcir = accumarray(ic, bin_table.DCIR, [], @mean);
    soc_centers = soc_bins(unique_bins(:,1));
    tbatt_centers = tbatt_bins(unique_bins(:,2));

    % scatter로 실제 데이터가 있는 bin만 표시
    figure('Position', [100, 100, 800, 600]);
    scatter(soc_centers, tbatt_centers, 200, mean_dcir, 's', 'filled');
    colormap(parula);
    colorbar;
    xlabel('SOC (%)');
    ylabel('Battery Temp. (°C)');
    title('DCIR in m\Omega (only bins with data)');
    set(gca, 'XTick', soc_bins, 'YTick', tbatt_bins, 'YDir', 'normal');
    grid on;

    % 셀 중앙에 값 표시
    for k = 1:length(mean_dcir)
        text(soc_centers(k), tbatt_centers(k), sprintf('%d', round(mean_dcir(k))), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
    end

    saveas(gcf, fullfile(saveDir, 'DCIR_SOC_Tbatt_heatmap_sparse.fig'));
    close(gcf);
end