function Psingle_and = trial_tot(D,Naverage,g,N)
% clear;clc
% Naverage=1; seed=22; D=0.005; g=5.484416576121011; N=100;  % 5.484416576121011   0.367466194073669

rng(1,"twister");
b=-0.5;

swtich_time_point=0.0;

Psingle_and = zeros(N, Naverage); 
stat_success_AND = zeros(N, Naverage);  % N × Naverage
stat_num = 0;

step = 0.01;  Tinterval = 1000;  T0 = 100;  Ttot = 7*Tinterval;
Nt = round(Ttot/step)+1;
t = 0:step:Ttot;

% xm = zeros(1, Naverage);
% xm_ = zeros(3, Naverage);
xm(1:3)=0;
k=0;
x = -1.0 * ones(N, Naverage);  % N × Naverage
%% 时间循环
for n  = 1:Nt
    if mod(n-1,round(Tinterval/step))==0 
            k=k+1;
            swtich_time_point = t(n);
            if k == 1 
                Iinput=1;
            elseif k == 2
                Iinput=0;
            elseif k == 3
                Iinput=-1;
            elseif k == 4 
                Iinput=0;
            elseif k == 5
                Iinput=1;
            elseif k == 6 
                Iinput=-1;
            elseif k == 7 
                Iinput=1;
            end
     end

    WhiteNoise=sqrt(-4.0*D*step.*log(rand(N,Naverage))).*cos(2.0*pi.*rand(N,Naverage));
    X_mean = mean(x,1);
    dx = x + step*(4*x-20*x.^3+b+Iinput+(g)*(X_mean - x)) + WhiteNoise;
    x = dx;
    
    % for k = 1:Naverage
    % xm_(1,k) = find_xm_only(-1, 1.5, 0.001, g, X_mean(k));
    % xm_(2,k) = find_xm_only(0, 1.5, 0.001, g, X_mean(k));
    % xm_(3,k) = find_xm_only(1, 1.5, 0.001, g, X_mean(k));
    % end
    % xm(Iinput==-1) = xm_(1, Iinput==-1);
    % xm(Iinput==0)  = xm_(2, Iinput==0);
    % xm(Iinput==1)  = xm_(3, Iinput==1);

    if t(n)-swtich_time_point>=T0
        stat_num = stat_num + 1;
      
        stat_success_AND(round(Iinput)==-1 & x<xm(1))=stat_success_AND(round(Iinput)==-1 & x<xm(1))+1;
        stat_success_AND(round(Iinput)==0 & x<xm(2))=stat_success_AND(round(Iinput)==0 & x<xm(2))+1;
        stat_success_AND(round(Iinput)==1 & x>xm(3))=stat_success_AND(round(Iinput)==1 & x>xm(3))+1;
    end
end
Psingle_and(stat_success_AND == stat_num) = 1; % 满足就赋值 1

% mean(Psingle_and,'all');