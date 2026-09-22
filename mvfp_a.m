%   图1 = 密度演化 / readout / 有效势 / x-P0 分叉
%   图2 = 冻结势近似下左右首次通过时间 MFPT vs P0
% 另含 MFPT vs D（固定 g,P0）与 PL(t) 左阱占比辅助图。
clear; clc; close all;

%% ==================== 主参数 ====================
% 局域双阱漂移: dx/dt = u*x - v*x^3 + ... + bias + noise
par.u      = 4;       
par.v      = 20;      
par.bias   = -0.5;    

par.D      = 0.03;   
par.sigma1 = 0.3;    
par.sigma2 = 0;       

opts.P0 = 0;         
opts.T = 100;        

opts.MFPT_P0 = 0;    
opts.MFPT_g = 0.3;    

opts.xMax = 1.0;     
opts.xGridN = 601;    
opts.initStd = 0.02;  
opts.odeRelTol = 1e-4;  
opts.odeAbsTol = 1e-8;
opts.thresholdMode = 'saddle'; 

opts.tSnapshots = [0 10 50 100];
opts.P0Scan = linspace(0, 1, 41);
opts.mfptP0Scan = linspace(0, 0.54, 28);
opts.potentialP0List = [0 0.25 0.5 0.75 1.0];
opts.DList = [0.005 0.01 0.02 0.05 0.1];
opts.sigma1List = [0 0.01 0.3 1 3 10];
opts.steadyFGrid = linspace(-6, 6, 901);
opts.steadyRootTol = 1e-8;

outDir = fullfile(pwd, 'outputs', 'student_minimal_mvfp');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
%% ==================== 计算 ====================
% 1) MVFP 密度轨迹 rho(x,t) + 稳态密度（自洽来自 M1,M2 平均场）
fprintf('Computing MVFP density trace...\n');
trace = solve_mvfp_density_trace(par, opts.P0, unique([0 opts.tSnapshots opts.T, 1000]), opts);
steadyDensity = selected_steady_density(par, opts);

% 2) 扫描 P0 得有限时间 readout P(T) + 稳态分支
fprintf('Computing finite-time readout curve...\n');
PByP0 = compute_readout_by_p0(par, opts);
steadyBranches = compute_steady_branches_for_param(par, opts);

% 3) 冻结平均场 F(P0) 下的有效势族 + x-P0 分叉
fprintf('Computing potential panels and x-P0 bifurcation...\n');
potD = compute_potential_vs_D(par, opts);
potP0 = compute_potential_vs_P0(par, opts);
potSigma1 = compute_potential_vs_sigma1(par, opts);
xBif = compute_x_bifurcation_vs_p0(par, opts);
% 导出 x-P0 分叉数据（Origin）
fid = fopen(fullfile(outDir, 'Figure_xP0_D0.03_g0.3.dat'), 'w');
fprintf(fid, 'P0\txL\txs\txR\txSingle\n');
for k = 1:length(xBif.P0)
    fprintf(fid, '%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n', ...
        xBif.P0(k), ...
        xBif.xL(k), ...
        xBif.xs(k), ...
        xBif.xR(k), ...
        xBif.xSingle(k));
end
fclose(fid);
fprintf('x-P0 bifurcation data saved:\n%s\n', ...
    fullfile(outDir, 'Figure_xP0_D0.03_g0.3.dat'));

% 4) 冻结势近似下 MFPT vs P0（方向不对称诊断量）
fprintf('Computing MFPT curves...\n');
mfpt = compute_mfpt_vs_p0(par, opts);
%% ==================== MFPT vs D（固定 g、P0） ====================
fprintf('Computing MFPT vs D at fixed g=%.5g, P0=%.2f...\n', ...
    opts.MFPT_g, opts.MFPT_P0);
parMFPT_D = par;
parMFPT_D.sigma1 = opts.MFPT_g;   
P0MFPT_D = opts.MFPT_P0;          
DminMFPT = 0.001;
DmaxMFPT = 0.3;
DScanMFPT = DminMFPT * ...        
    ((DmaxMFPT / DminMFPT).^(1/95)).^(0:95);
mfptD = compute_mfpt_vs_D( ...
    parMFPT_D, P0MFPT_D, DScanMFPT, opts);
%% -------------------- Plot MFPT vs D --------------------
figMFPT_D = figure('Color','w', ...
    'Name','MFPT vs D at fixed g and P0', ...
    'Units','pixels', ...
    'Position',[150 150 900 600]);
semilogx(mfptD.D, mfptD.log10TauLR, ...
    'bo-', 'LineWidth',1.3, 'MarkerSize',5);
hold on;

semilogx(mfptD.D, mfptD.log10TauRL, ...
    'ro-', 'LineWidth',1.3, 'MarkerSize',5);
xlabel('D');
ylabel('log_{10}(\tau)');
title(sprintf('MFPT vs D, g=%.4g, P_0=%.2f', ...
    parMFPT_D.sigma1, P0MFPT_D));
legend('\tau_{L\rightarrow R}', ...
       '\tau_{R\rightarrow L}', ...
       'Location','best');
grid on;
box on;
%% -------------------- Export MFPT-D data --------------------
fid = fopen(fullfile(outDir, 'MFPT_vs_D_g0.3.dat'), 'w');
fprintf(fid, ...
    'D\tlogD\tlog10TauLR\tlog10TauRL\tlog10TauRatio\n');
for k = 1:length(mfptD.D)
    fprintf(fid, ...
        '%.12g\t%.12f\t%.12g\t%.12g\t%.12g\n', ...
        mfptD.D(k), ...
        log10(mfptD.D(k)), ...
        mfptD.log10TauLR(k), ...
        mfptD.log10TauRL(k), ...
        mfptD.log10TauRatio(k));
end
fclose(fid);
fprintf('MFPT-D data saved to:\n%s\n', ...
    fullfile(outDir, 'MFPT_vs_D_g0.3.dat'));
%% ==================== 图1: 核心 MVFP 机制 ====================
fig1 = figure('Color', 'w', 'Name', 'Minimal MVFP core', ...
    'Units', 'pixels', 'Position', [80 80 1350 950]);
tiledlayout(fig1, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot_density_panel(trace, steadyDensity, opts);                                   % 密度演化

nexttile;
plot_potential_family(potD.x, potD.U, potD.labels, 'Potential vs D');             % 势 vs D

nexttile;
plot_readout_panel(PByP0, steadyBranches, opts);                                  % readout P(T|P0)

nexttile;
plot_potential_family(potP0.x, potP0.U, potP0.labels, 'Potential vs P_0');        % 势 vs P0

nexttile;
plot_potential_family(potSigma1.x, potSigma1.U, potSigma1.labels, 'Potential vs \sigma_1'); % 势 vs sigma1

nexttile;
plot_x_p0_bifurcation(xBif);                                                      % x-P0 分叉

%% 导出 Potential vs sigma1 每条减自身最小值
fid2 = fopen(fullfile(outDir,'Potential_vs_sigma3.dat'),'w');
fprintf(fid2,'# x\t');
fprintf(fid2,'%s\t',potSigma1.labels{:}); fprintf(fid2,'\n');
U_shift = potSigma1.U - min(potSigma1.U, [], 1);
data = [potSigma1.x, U_shift];
for i=1:size(data,1)
    fprintf(fid2,'%.6f\t',data(i,:)); fprintf(fid2,'\n');
end
fclose(fid2);

save_figure(fig1, fullfile(outDir, 'student_fig_mvfp_core_self_consistent_switch.png'));
%% ==================== 图2: MFPT vs P0 ====================
%% 导出 MFPT 数据
fid = fopen(fullfile(outDir, 'Figure_MFPT_D0.03_g0.3.dat'), 'w');
fprintf(fid, 'P0\tlog10TauLR\tlog10TauRL\tlog10TauRatio\n');
log10Ratio = (mfpt.logTauRL - mfpt.logTauLR) ./ log(10);
for k = 1:length(mfpt.P0)
    fprintf(fid, '%.12g\t%.12g\t%.12g\t%.12g\n', ...
        mfpt.P0(k), ...
        mfpt.log10TauLR(k), ...
        mfpt.log10TauRL(k), ...
        log10Ratio(k));
end
fclose(fid);
fprintf('Figure 2 MFPT data saved:\n%s\n', ...
    fullfile(outDir, 'Figure_MFPT_D0.03_g0.3.dat'));

fig2 = figure('Color', 'w', 'Name', 'Minimal MVFP MFPT', ...
    'Units', 'pixels', 'Position', [120 120 900 620]);
tiledlayout(fig2, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
idx = isfinite(mfpt.log10TauLR);
plot(mfpt.P0(idx), mfpt.log10TauLR(idx), 'b-o', 'LineWidth', 1.3);
xlabel('P_0'); ylabel('log_{10}\tau_{L\rightarrow R}');
title('First passage from left well'); grid on;

nexttile;
idx = isfinite(mfpt.log10TauRL);
plot(mfpt.P0(idx), mfpt.log10TauRL(idx), 'r-o', 'LineWidth', 1.3);
xlabel('P_0'); ylabel('log_{10}\tau_{R\rightarrow L}');
title('First passage from right well'); grid on;

nexttile;
log10Ratio = (mfpt.logTauRL - mfpt.logTauLR) ./ log(10);
idx = isfinite(log10Ratio);
plot(mfpt.P0(idx), log10Ratio(idx), 'k-o', 'LineWidth', 1.3);
xlabel('P_0'); ylabel('log_{10}(\tau_{R\rightarrow L}/\tau_{L\rightarrow R})');
title(sprintf('Directional asymmetry, D=%.3g', par.D)); grid on;

save_figure(fig2, fullfile(outDir, 'student_fig_mvfp_mfpt_vs_p0.png'));

%% ==================== PL(t) 曲线: x<0 左阱粒子占比 (0~1000s) ====================
fprintf('Generating PL(t) left-well fraction figure (0~1000s)...\n');
gridPL = make_grid(opts);
xGrid = gridPL.x;
xLeftMask = xGrid < 0;   % x<0 判定为左阱

maskTimeWin = (trace.tFull >= 0) & (trace.tFull <= 1000);
tPlotPL = trace.tFull(maskTimeWin);
nPL = length(tPlotPL);
PLplot = zeros(nPL, 1);
% 逐时刻积分左阱占比
for kt = 1:nPL
    tNow = tPlotPL(kt);
    [~, idxT] = min(abs(trace.tFull - tNow));
    rhoNow = trace.rhoAll(:, idxT);
    massLeft = trapz(xGrid(xLeftMask), rhoNow(xLeftMask));
    PLplot(kt) = min(max(massLeft, 0), 1);
end
dt_fixed = 1;
t_ori = tPlotPL;
PL_ori = PLplot;
tPlotPL = linspace(0, 1000, round(1000/dt_fixed)+1);
PLplot = interp1(t_ori, PL_ori, tPlotPL, 'linear','extrap');

figPL = figure('Color', 'w', 'Name', 'PL(t) Left-well fraction x<0',...
    'Units','pixels','Position',[150 150 850 520]);
plot(tPlotPL, PLplot, 'k-', 'LineWidth',1.5);
xlabel('t');
ylabel('P_L(t) (Fraction of particles x < 0)');
title('Left-well population ratio over time (0 \le t \le 1000)','Interpreter','none');
grid on;
ylim([0, 1.02]);
xlim([0, 1000]);
save_figure(figPL, fullfile(outDir, 'fig_PL_time_evolution.png'));

fprintf('Done. Figures saved under:\n  %s\n', outDir);

%% 绘图辅助函数
function plot_density_panel(trace, steadyDensity, opts)
% 密度面板：rho(x,t) 各时刻与稳态 rho_ss 对比
hold on;
for k = 1:numel(trace.t)
    plot(trace.x, trace.rho(:, k), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('t=%g', trace.t(k)));
end
if ~isempty(steadyDensity.rho)
    plot(steadyDensity.x, steadyDensity.rho, 'k-', ...
        'LineWidth', 2.0, 'DisplayName', '\rho_{ss}');
end
xline(trace.threshold(end), '--', 'Color', 0.55 .* [1 1 1], ...
    'DisplayName', 'final threshold');
hold off;
xlabel('x'); ylabel('\rho(x,t)');
title(sprintf('Density, P_0=%g', opts.P0));
grid on; legend('Location', 'best');
end

function plot_readout_panel(PByP0, steadyBranches, opts)
% readout 面板：P(T|P0) 是否存在有限时间集体切换
plot(PByP0.P0, PByP0.PT, 'o-', 'LineWidth', 1.4, ...
    'DisplayName', 'finite time');
hold on;
stable = logical(steadyBranches.stable);
if any(stable)
    for k = find(stable)
        yline(steadyBranches.Pss(k), 'k-', 'LineWidth', 1.2, ...
            'DisplayName', 'steady stable');
    end
end
plot([0 1], [0 1], '--', 'Color', 0.65 .* [1 1 1], ...
    'DisplayName', 'P(T)=P_0');
hold off;
xlabel('initial P_0'); ylabel(sprintf('P(t=%g)', opts.T));
title('Finite-time readout'); ylim([0 1.05]);
grid on; legend('Location', 'best');
end

function plot_potential_family(x, U, labels, titleText)
% 势能族：每条减自身最小值，比较势垒形状
hold on;
for k = 1:size(U, 2)
    Uk = U(:, k);
    plot(x, Uk - min(Uk), 'LineWidth', 1.3, 'DisplayName', labels{k});
end
hold off;
xlabel('x'); ylabel('\Delta U_{eff}');
title(titleText); xlim([-1 1]); ylim([-0.1 1.0]);
grid on; legend('Location', 'best');
end

function plot_x_p0_bifurcation(bif)
% x-P0 分叉：冻结势极值位置随 P0 变化（双阱/单阱）
hold on;
plot(bif.P0, bif.xL, '-', 'LineWidth', 1.4, 'DisplayName', 'left minimum');
plot(bif.P0, bif.xs, '--', 'LineWidth', 1.2, 'DisplayName', 'saddle');
plot(bif.P0, bif.xR, '-', 'LineWidth', 1.4, 'DisplayName', 'right minimum');
plot(bif.P0, bif.xSingle, 'k-', 'LineWidth', 1.4, 'DisplayName', 'single minimum');
hold off;
xlabel('P_0'); ylabel('x'); title('x-P_0 bifurcation');
xlim([min(bif.P0), max(bif.P0)]); ylim([-0.5 0.4]);
grid on; legend('Location', 'best');
end

function save_figure(fig, path)
if exist('exportgraphics', 'file')
    exportgraphics(fig, path, 'Resolution', 200);
else
    saveas(fig, path);
end
end

%% 核心 MVFP 密度求解
function trace = solve_mvfp_density_trace(par, P0, times, opts)
% 求解一维 MVFP: dt rho = -dx[b(x,rho) rho] + D dxx rho
% 平均场反馈来自矩 M1=mean(x), M2=mean(x^2)
times = unique(times(:).');
grid = make_grid(opts);
x = grid.x;
rho0 = initial_density(x, par, P0, opts.initStd);
tSolve = unique([0 times]);
if max(tSolve) <= 0
    tOut = 0;
    y = rho0.';
else
    odeOpts = odeset('RelTol', opts.odeRelTol, 'AbsTol', opts.odeAbsTol, ...
        'NonNegative', 1:numel(x));
    rhs = @(~, rho) mvfp_rhs(grid, rho, par);
    [tOut, y] = ode15s(rhs, tSolve, rho0, odeOpts);
end

fullPL = zeros(size(tOut));
fullRhoAll = zeros(numel(x), numel(tOut));
thresholdAll = zeros(size(tOut));
for kt = 1:numel(tOut)
    rho = y(kt, :).';
    rho = normalize_density(x, rho);
    thresholdAll(kt) = density_threshold(x, rho, par, opts.thresholdMode);
    leftMask = x < 0;
    fullPL(kt) = min(max(trapz(x(leftMask), rho(leftMask)), 0), 1);
    fullRhoAll(:,kt) = rho;
end

rhoOut = nan(numel(x), numel(times));
P = nan(numel(times), 1);
threshold = nan(numel(times), 1);
rawMass = nan(numel(times), 1);
for k = 1:numel(times)
    [~, idx] = min(abs(tOut - times(k)));
    rho = y(idx, :).';
    rawMass(k) = trapz(x, rho);
    rho = normalize_density(x, rho);
    threshold(k) = density_threshold(x, rho, par, opts.thresholdMode);
    left = x < 0;
    P(k) = min(max(trapz(x(left), rho(left)), 0), 1);
    rhoOut(:, k) = rho;
end
trace = struct('x', x, 't', times, 'rho', rhoOut, ...
    'P', P, 'threshold', threshold, 'rawMass', rawMass,...
    'tFull', tOut, 'PLfull', fullPL, 'rhoAll', fullRhoAll);
end

function drho = mvfp_rhs(grid, rho, par)
% MVFP 右端项：通量 J = b*rho - D*drho/dx，零通量边界保总质量
x = grid.x;
rho = rho(:);
mass = grid.w.' * rho;
if mass > 0 && isfinite(mass)
    M1 = (grid.wx.' * rho) ./ mass;
    M2 = (grid.wx2.' * rho) ./ mass;
else
    M1 = grid.wx.' * rho;
    M2 = grid.wx2.' * rho;
end
b = (par.u - par.sigma1) .* x - (par.v + par.sigma2) .* grid.x3 ...
    + par.sigma1 .* M1 + par.sigma2 .* M2 .* M1 + par.bias;
bFace = 0.5 .* (b(1:end-1) + b(2:end));
rhoLeft = rho(1:end-1);
rhoRight = rho(2:end);
rhoUpwind = rhoLeft;
rhoUpwind(bFace < 0) = rhoRight(bFace < 0);
Jadv = bFace .* rhoUpwind;
Jdiff = -par.D .* (rhoRight - rhoLeft) ./ grid.dx;
J = Jadv + Jdiff;
Jfull = [0; J; 0];
drho = -(Jfull(2:end) - Jfull(1:end-1)) ./ grid.dx;
end

function grid = make_grid(opts)
% 一维均匀网格 + 梯形积分权重
x = linspace(-opts.xMax, opts.xMax, opts.xGridN).';
dx = x(2) - x(1);
w = ones(size(x)) .* dx;
w([1 end]) = dx ./ 2;
grid = struct('x', x, 'dx', dx, 'w', w, ...
    'wx', w .* x, 'wx2', w .* x.^2, 'x3', x.^3);
end

function rho = initial_density(x, par, P0, initStd)
% 初始双高斯包：P0 控制左阱权重
info0 = frozenP_initial_info(P0, par);
[xL, xR] = initial_representatives(info0, par);
rho = P0 .* gaussian_pdf(x, xL, initStd) ...
    + (1 - P0) .* gaussian_pdf(x, xR, initStd);
rho = normalize_density(x, rho);
end

function info = frozenP_initial_info(P0, par)
% 冻结平均场 F(P0) 下的左右极小值
info = struct('valid', false, 'xLstar', nan, 'xRstar', nan);
F = initial_field_guess(P0, par);
[ok, xL, ~, xR] = extrema_from_F(F, par);
if ok
    info.valid = true;
    info.xLstar = xL;
    info.xRstar = xR;
end
end

function [xLinit, xRinit] = initial_representatives(info, par)
% 初始高斯中心：双阱极小值，否则解析近似
if info.valid
    xLinit = info.xLstar;
    xRinit = info.xRstar;
    return;
end
[ok, xL, ~, xR] = extrema_from_F(par.bias, par);
if ok
    xLinit = xL;
    xRinit = xR;
    return;
end
A = par.u - par.sigma1;
B = par.v + par.sigma2;
if A > 0 && B > 0
    x0 = sqrt(max(A ./ B, eps));
else
    x0 = 1;
end
xLinit = -x0;
xRinit = x0;
end

function y = gaussian_pdf(x, mu, sigma)
y = exp(-0.5 .* ((x - mu) ./ sigma).^2) ./ (sqrt(2 .* pi) .* sigma);
end

function rho = normalize_density(x, rho)
% 归一化并去掉负值
rho = max(rho(:), 0);
mass = trapz(x, rho);
if mass > 0 && isfinite(mass)
    rho = rho ./ mass;
else
    rho(:) = 0;
    [~, idx] = min(abs(x));
    rho(idx) = 1 ./ max(eps, trapz(x, double((1:numel(x)).' == idx)));
end
end

function threshold = density_threshold(x, rho, par, mode)
% readout 分界线：瞬时有效势鞍点，否则回退 x=0
if strcmpi(mode, 'saddle')
    M1 = trapz(x, x .* rho);
    M2 = trapz(x, x.^2 .* rho);
    F = par.sigma1 .* M1 + par.sigma2 .* M2 .* M1 + par.bias;
    [ok, ~, xs, ~] = extrema_from_F(F, par);
    if ok
        threshold = xs;
        return;
    end
end
threshold = 0;
end

%% Readout / 势能 / 分叉面板
function PByP0 = compute_readout_by_p0(par, opts)
% 扫描 P0，每个 P0 完整求解 MVFP，T 时刻读左侧概率
P0Vals = opts.P0Scan(:).';
PT = nan(size(P0Vals));
for k = 1:numel(P0Vals)
    tr = solve_mvfp_density_trace(par, P0Vals(k), opts.T, opts);
    PT(k) = tr.P(end);
end
PByP0 = struct('P0', P0Vals, 'PT', PT);
end

function pot = compute_potential_vs_D(par, opts)
% 冻结有效势 vs D
x = linspace(-opts.xMax, opts.xMax, opts.xGridN).';
DVals = opts.DList(:).';
U = nan(numel(x), numel(DVals));
labels = cell(size(DVals));
for k = 1:numel(DVals)
    p = par;
    p.D = DVals(k);
    F = initial_field_guess(opts.P0, p);
    U(:, k) = effective_potential_from_F(x, p, F);
    labels{k} = sprintf('D=%g', DVals(k));
end
pot = struct('x', x, 'U', U, 'labels', {labels});
end

function pot = compute_potential_vs_P0(par, opts)
% 冻结有效势 vs P0（P0 改变初始平均场倾斜）
x = linspace(-opts.xMax, opts.xMax, opts.xGridN).';
P0Vals = opts.potentialP0List(:).';
U = nan(numel(x), numel(P0Vals));
labels = cell(size(P0Vals));
for k = 1:numel(P0Vals)
    F = initial_field_guess(P0Vals(k), par);
    U(:, k) = effective_potential_from_F(x, par, F);
    labels{k} = sprintf('P_0=%g', P0Vals(k));
end
pot = struct('x', x, 'U', U, 'labels', {labels});
end

function pot = compute_potential_vs_sigma1(par, opts)
% 冻结有效势 vs sigma1
x = linspace(-opts.xMax, opts.xMax, opts.xGridN).';
sVals = opts.sigma1List(:).';
U = nan(numel(x), numel(sVals));
labels = cell(size(sVals));
for k = 1:numel(sVals)
    p = par;
    p.sigma1 = sVals(k);
    F = initial_field_guess(opts.P0, p);
    U(:, k) = effective_potential_from_F(x, p, F);
    labels{k} = sprintf('\\sigma_1=%g', sVals(k));
end
pot = struct('x', x, 'U', U, 'labels', {labels});
end

function bif = compute_x_bifurcation_vs_p0(par, opts)
% 冻结势极值位置随 P0 变化（x-P0 分叉）
P0Vals = opts.P0Scan(:).';
xL = nan(size(P0Vals));
xs = nan(size(P0Vals));
xR = nan(size(P0Vals));
xSingle = nan(size(P0Vals));
for k = 1:numel(P0Vals)
    F = initial_field_guess(P0Vals(k), par);
    extrema = effective_potential_extrema_from_F(F, par);
    xL(k) = extrema.xL;
    xs(k) = extrema.xs;
    xR(k) = extrema.xR;
    xSingle(k) = extrema.xSingle;
end
bif = struct('P0', P0Vals, 'xL', xL, 'xs', xs, ...
    'xR', xR, 'xSingle', xSingle);
end

%% 稳态密度与自洽分支
function steady = selected_steady_density(par, opts)
% 从自洽稳态中选一个稳定分支作密度对照（优先左阱概率高者）
x = linspace(-opts.xMax, opts.xMax, opts.xGridN).';
branches = compute_steady_branches_for_param(par, opts);
if isempty(branches.F)
    steady = struct('x', x, 'rho', []);
    return;
end
stable = logical(branches.stable);
idx = find(stable);
if isempty(idx)
    idx = 1:numel(branches.F);
end
[~, pickLocal] = max(branches.Pss(idx));
pick = idx(pickLocal);
[rho, ~, ~, ~] = steady_density_from_F(branches.F(pick), par, x);
steady = struct('x', x, 'rho', rho);
end

function branches = compute_steady_branches_for_param(par, opts)
% 稳态自洽方程 F = sigma1*M1 + sigma2*M2*M1 + bias 求根

x = linspace(-opts.xMax, opts.xMax, opts.xGridN).';
FGrid = opts.steadyFGrid(:).';
residualGrid = nan(size(FGrid));
for k = 1:numel(FGrid)
    residualGrid(k) = steady_residual(FGrid(k), par, x);
end
tol = opts.steadyRootTol;
rootCandidates = [];
nearMask = isfinite(residualGrid) & abs(residualGrid) <= tol;
rootCandidates = [rootCandidates, FGrid(nearMask)]; 
for k = 1:(numel(FGrid) - 1)
    r1 = residualGrid(k);
    r2 = residualGrid(k + 1);
    if isfinite(r1) && isfinite(r2) && r1 .* r2 < 0
        try
            rootCandidates(end + 1) = fzero(@(F) steady_residual(F, par, x), ...
                [FGrid(k), FGrid(k + 1)]); 
        catch
        end
    end
end
rootsF = unique_roots(rootCandidates, max(10 .* tol, 1e-7));
F = nan(1, numel(rootsF));
m1 = F; m2 = F; Pss = F; stable = F; residualAtRoot = F;
for k = 1:numel(rootsF)
    F(k) = rootsF(k);
    [residualAtRoot(k), m1(k), m2(k), Pss(k)] = steady_residual(F(k), par, x);
    stable(k) = steady_branch_stability(F(k), par, x);
end
branches = struct('F', F, 'm1', m1, 'm2', m2, ...
    'Pss', Pss, 'stable', stable, 'residualAtRoot', residualAtRoot);
end

function [residual, m1, m2, Pss] = steady_residual(F, par, x)
% 自洽残差：F 与平均场反馈之差
[~, m1, m2, Pss] = steady_density_from_F(F, par, x);
feedback = par.sigma1 .* m1 + par.sigma2 .* m2 .* m1 + par.bias;
residual = F - feedback;
end

function [rho, m1, m2, Pss] = steady_density_from_F(F, par, x)
% 冻结 F 的 Boltzmann 稳态密度 rho_ss ~ exp(-U_eff/D)
U = effective_potential_from_F(x, par, F);
logWeight = -U ./ par.D;
shift = max(logWeight);
w = exp(logWeight - shift);
Z = trapz(x, w);
rho = w ./ Z;
m1 = trapz(x, x .* rho);
m2 = trapz(x, x.^2 .* rho);
[ok, ~, xs, ~] = extrema_from_F(F, par);
if ok
    threshold = xs;
else
    threshold = 0;
end
Pss = integrate_density_left_of_threshold(x, rho, threshold);
end

function P = integrate_density_left_of_threshold(x, rho, threshold)
% 阈值左侧概率（限制到 [0,1]）
left = x <= threshold;
P = min(max(trapz(x(left), rho(left)), 0), 1);
end

function stable = steady_branch_stability(F, par, x)
% 分支稳定性：固定点映射斜率 dG/dF < 1
dF = max(1e-5, 1e-4 * max(1, abs(F)));
[~, m1p, m2p] = steady_residual(F + dF, par, x);
[~, m1m, m2m] = steady_residual(F - dF, par, x);
Gp = par.sigma1 .* m1p + par.sigma2 .* m2p .* m1p + par.bias;
Gm = par.sigma1 .* m1m + par.sigma2 .* m2m .* m1m + par.bias;
dGdF = (Gp - Gm) ./ (2 .* dF);
stable = double(isfinite(dGdF) && dGdF < 1);
end

function rootsOut = unique_roots(values, tol)
% 合并重复根
values = sort(values(isfinite(values)));
rootsOut = [];
for k = 1:numel(values)
    if isempty(rootsOut) || abs(values(k) - rootsOut(end)) > tol
        rootsOut(end + 1) = values(k); %#ok<AGROW>
    end
end
end

%% 有效势与极值
function F = initial_field_guess(P, par)
% 初始平均场估计: M1=P*xL+(1-P)*xR, M2=P*xL^2+(1-P)*xR^2
A = par.u - par.sigma1;
B = par.v + par.sigma2;
x0 = sqrt(max(A ./ B, eps));
xL = -x0;
xR = x0;
[valid, x1, ~, x2] = extrema_from_F(par.bias, par);
if valid
    xL = x1;
    xR = x2;
end
M1 = P .* xL + (1 - P) .* xR;
M2 = P .* xL.^2 + (1 - P) .* xR.^2;
F = par.sigma1 .* M1 + par.sigma2 .* M2 .* M1 + par.bias;
end

function U = effective_potential_from_F(x, par, F)
% 有效势: U = -(u-s1)x^2/2 + (v+s2)x^4/4 - F*x（drift = -dU/dx）
A = par.u - par.sigma1;
B = par.v + par.sigma2;
U = -0.5 .* A .* x.^2 + 0.25 .* B .* x.^4 - F .* x;
end

function extrema = effective_potential_extrema_from_F(F, par)
% 冻结势极值：双阱(左极小/鞍点/右极小) 或 单阱(single)
[ok, xL, xs, xR] = extrema_from_F(F, par);
extrema = struct('xL', nan, 'xs', nan, 'xR', nan, 'xSingle', nan);
if ok
    extrema.xL = xL;
    extrema.xs = xs;
    extrema.xR = xR;
    return;
end
A = par.u - par.sigma1;
B = par.v + par.sigma2;
r = roots([B, 0, -A, -F]);
r = sort(real(r(abs(imag(r)) < 1e-10)));
if numel(r) == 1
    extrema.xSingle = r(1);
end
end

function [ok, xL, xs, xR] = extrema_from_F(F, par)
% 极值满足 B*x^3 - A*x - F = 0；曲率 +,-,+ 对应 左极小/鞍点/右极小
A = par.u - par.sigma1;
B = par.v + par.sigma2;
ok = false; xL = nan; xs = nan; xR = nan;
if ~(isfinite(A) && isfinite(B) && isfinite(F)) || B <= 0
    return;
end
r = roots([B, 0, -A, -F]);
r = sort(real(r(abs(imag(r)) < 1e-10)));
if numel(r) < 3
    return;
end
curv = -A + 3 .* B .* r.^2;
if curv(1) > 0 && curv(2) < 0 && curv(3) > 0
    ok = true;
    xL = r(1);
    xs = r(2);
    xR = r(3);
end
end

%% MFPT：backward Fokker-Planck log 积分
function mfpt = compute_mfpt_vs_p0(par, opts)
% 每个 P0 构造冻结势，算左→右、右→左平均首次通过时间
xMin = -opts.xMax;
xMax = opts.xMax;
domainGridN = max(101, opts.xGridN);
P0Vals = opts.mfptP0Scan(:).';
logTauLR = nan(size(P0Vals));
logTauRL = nan(size(P0Vals));
for k = 1:numel(P0Vals)
    F = initial_field_guess(P0Vals(k), par);
    [ok, ~, xs, ~] = extrema_from_F(F, par);
    if ~ok
        continue;
    end
    [~, logTauLR(k)] = mfpt_logquad_average_to_saddle( ...
        xMin, xs, par, F, domainGridN, 'right');
    [~, logTauRL(k)] = mfpt_logquad_average_to_saddle( ...
        xs, xMax, par, F, domainGridN, 'left');
end
mfpt = struct('P0', P0Vals, 'logTauLR', logTauLR, ...
    'logTauRL', logTauRL, 'log10TauLR', logTauLR ./ log(10), ...
    'log10TauRL', logTauRL ./ log(10));
end

function mfpt = compute_mfpt_vs_D(par, P0, DVals, opts)
% 固定 g=P0，扫描 D 的左右方向 MFPT
xMin = -opts.xMax;          xMax = opts.xMax;
domainGridN = max(101, opts.xGridN);
DVals = DVals(:).';
logTauLR = nan(size(DVals));   logTauRL = nan(size(DVals));
% 冻结平均场 F 与 D 无关
F = initial_field_guess(P0, par);
[ok, ~, xs, ~] = extrema_from_F(F, par);

if ~ok
    warning(['当前 g=%.5g, P0=%.5g 下冻结有效势没有双阱结构，' ...
             '无法计算 MFPT。'], ...
             par.sigma1, P0);
    mfpt = struct( ...
        'D', DVals, ...
        'logTauLR', logTauLR, ...
        'logTauRL', logTauRL, ...
        'log10TauLR', logTauLR ./ log(10), ...
        'log10TauRL', logTauRL ./ log(10), ...
        'log10TauRatio', nan(size(DVals)));
    return;
end
for k = 1:numel(DVals)
    p = par;
    p.D = DVals(k);
    [~, logTauLR(k)] = ...
        mfpt_logquad_average_to_saddle( ...
        xMin, xs, p, F, domainGridN, 'right');
    [~, logTauRL(k)] = ...
        mfpt_logquad_average_to_saddle( ...
        xs, xMax, p, F, domainGridN, 'left');
end
log10TauLR = logTauLR ./ log(10);
log10TauRL = logTauRL ./ log(10);
log10TauRatio = log10TauRL - log10TauLR;
mfpt = struct( 'D', DVals,     'logTauLR', logTauLR, ...
    'logTauRL', logTauRL,     'log10TauLR', log10TauLR, ...
    'log10TauRL', log10TauRL,   'log10TauRatio', log10TauRatio);
end

function [tauAvg, logTauAvg] = mfpt_logquad_average_to_saddle( ...
    xLeft, xRight, par, F, domainGridN, absorbSide)
% 一维 backward FP 积分公式求 MFPT，全程 log 域防溢出
tauAvg = nan; logTauAvg = nan;
if ~(isfinite(xLeft) && isfinite(xRight) && xRight > xLeft)
    return;
end
x = linspace(xLeft, xRight, max(5, round(domainGridN))).';
U = effective_potential_from_F(x, par, F);
logWeight = -U ./ par.D;
logZ = log_trapz_exp(x, logWeight);
if strcmpi(absorbSide, 'right')
    logInner = log_cumtrapz_exp(x, logWeight);
    logOuter = U ./ par.D + logInner;
    logTauProfile = log_cumtrapz_exp_reverse(x, logOuter) - log(par.D);
elseif strcmpi(absorbSide, 'left')
    logInner = log_cumtrapz_exp_reverse(x, logWeight);
    logOuter = U ./ par.D + logInner;
    logTauProfile = log_cumtrapz_exp(x, logOuter) - log(par.D);
else
    return;
end
logNumerator = log_trapz_exp(x, logTauProfile + logWeight);
logTauAvg = logNumerator - logZ;
tauAvg = exp_if_representable(logTauAvg);
end

function logVal = log_trapz_exp(x, logf)
% log(int exp(logf) dx)：log-sum-exp 梯形积分
x = x(:); logf = logf(:);
terms = -inf(numel(x) - 1, 1);
for i = 1:(numel(x) - 1)
    dx = x(i + 1) - x(i);
    terms(i) = log(dx ./ 2) + log_add_exp(logf(i), logf(i + 1));
end
logVal = -inf;
for i = 1:numel(terms)
    logVal = log_add_exp(logVal, terms(i));
end
end

function logCum = log_cumtrapz_exp(x, logf)
% 左→右累计 log 梯形积分
x = x(:); logf = logf(:);
logCum = -inf(numel(x), 1);
acc = -inf;
for i = 1:(numel(x) - 1)
    dx = x(i + 1) - x(i);
    term = log(dx ./ 2) + log_add_exp(logf(i), logf(i + 1));
    acc = log_add_exp(acc, term);
    logCum(i + 1) = acc;
end
end

function logCum = log_cumtrapz_exp_reverse(x, logf)
% 右→左累计 log 梯形积分（反向吸收边界）
x = x(:); logf = logf(:);
logCum = -inf(numel(x), 1);
acc = -inf;
for i = numel(x):-1:2
    dx = x(i) - x(i - 1);
    term = log(dx ./ 2) + log_add_exp(logf(i), logf(i - 1));
    acc = log_add_exp(acc, term);
    logCum(i - 1) = acc;
end
end

function c = log_add_exp(a, b)
% 稳定计算 log(exp(a)+exp(b))
if ~isfinite(a)
    c = b; return;
end
if ~isfinite(b)
    c = a; return;
end
m = max(a, b);
c = m + log(exp(a - m) + exp(b - m));
end

function y = exp_if_representable(logY)
% logY 未超 double 上限才转回 exp
if logY < log(realmax)
    y = exp(logY);
else
    y = inf;
end
end
