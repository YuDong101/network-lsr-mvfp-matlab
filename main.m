clear;clc
rng(1);
%% 参数设置
Naverage=1;
Dmin=0.001;Dmax=0.3;
D = Dmin * ((Dmax/Dmin).^(1/95)) .^(0:95); 
gmin=0.01;gmax=10;
g_=gmin * ((gmax/gmin).^(1/95)) .^(0:95);
% g_= 5.25;
% g_=1:0.05:2;
N = 1000;  
Pand_d_and = zeros(length(g_), length(D), N);
%% 主计算循环
tic
for j=1:length(g_)
    g=g_(j);
    parfor d=1:length(D)
         P_single_and = trial_tot(D(d),Naverage,g,N);
        % Pand_d_and(j,d,:) = P_single_and(j,d,:);
        Pand_d_and(j,d,:) = P_single_and;
    end
    fprintf("j/%d 耗时 %.4f \n",j,toc)
   
     % save Al_N=1000_24_AND.mat
end

Pand_and = mean(Pand_d_and, 3);
  save Al_N=1000_96_AND.mat
%%
 % figure(32)
 % plot(log10(D), Pand_and)
%% 保存结果到文本文件
fp7 = fopen('N=1000_DgP_96.dat','w');
fprintf(fp7, 'log10(D) g log10(g) AND\n');
for jj=1:length(g_)
    for ii=1:length(D)         
       fprintf(fp7,'%f %f %f %f\n',log10(D(ii)),g_(jj),log10(g_(jj)),Pand_and(jj,ii));
    end
end
fclose(fp7);
%% 可视化AND门结果 - 热力图（imagesc）
figure(11)
set(gcf, 'Position', [100, 100, 500, 400]); 
imagesc(log10(D), log10(g_), Pand_and);   % Y 轴使用 log10(g_)
set(gca, 'YDir', 'normal');  % Y轴正向朝上
xlabel('log_{10}(D) - 噪声强度对数');
ylabel('log_{10}(g) - 全局耦合强度对数');
title(['全耦合AND门模型 (N=', num2str(N), ' 个逻辑门) - 成功概率热力图']);
colorbar;
colormap(jet);
clim([0 1]);  % 固定颜色范围 0-1
