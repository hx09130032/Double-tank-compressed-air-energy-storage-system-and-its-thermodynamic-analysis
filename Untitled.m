%2020.12.19
%% （单级压缩机+等温循环压缩）/（单级膨胀机+多变循环膨胀）
clc;
clear all;
close all;
for i=1:11
    for j=1:11
%% 以发电规模定储气罐容积大小
        % 储气罐压力 4-6MPa
        % 以 4MPa 压力和 Th 温度进行1MW规模释能1h(3.6*10^9 J)
        
k=1.4;                                       % 绝热指数
Rg=287;                                      % 气体常数
Tamb=273.15+20;                              % 环境温度 K
P_amb=101325;                                % 环境压力 Pa
Th=273.15+300;                               % 太阳能加热后气体温度 K
n_T=0.85;                                    % 透平效率
P_8=5500000+100000*(j-1);                                 % 储气罐最高压力 Pa
P_6=3500000+100000*(i-1);                                 % 储气罐最低压力 Pa
sss=sqrt(P_6/P_amb);
T_OUT=Th*(1-n_T*(1-(1/sss)^((k-1)/k)));     % 透平出口温度 K
w=k*Rg*Th/(k-1)*n_T*(1-(1/sss)^((k-1)/k));  % 单位质量单级膨胀做功 J
M=(3.6*10^9)/(2*w);                          % 气体总质量 kg
m=M/3600;                                    % 气体质量流量 kg/s
         % P_8*V-P_6*V=M*Rg*Tamb
V=M*Rg*Tamb/(P_8-P_6);                       % 储气罐容积 m^3
%% 储能过程
Vs=5;                                                   % 循环压缩装置单次压气量 m^3
V2=V;                                                   % 储气罐容积 m^3
n_C=0.85;                                               % 压缩机效率
          % 一级压缩         
P_xun_in=400000;        

epsilon=(P_xun_in/P_amb);                               % 压比
T_C1_out=Tamb+[(epsilon)^((k-1)/k)-1]*Tamb/n_C;         % 第一级压缩机排气温度 K

T_xun_in=Tamb;                                          % 循环压缩装置进气温度 K

d=P_xun_in*Tamb*Vs/T_xun_in/V2;                         % 循环压缩单次增压量 Pa
P_8_end=(fix((P_8-P_6)/d)+1)*d+P_6;                     % 储气罐真实终了压力 Pa
n=(P_8_end-P_6)/d;                                      % 压缩次数
w_C1=(1/n_C)*Rg*Tamb*(k/(k-1))*[epsilon^((k-1)/k)-1];   % 单位质量一级压缩耗功 J
M_TOTAL=P_8_end*V2/Rg/Tamb-P_6*V2/Rg/Tamb;              % 压缩过程气体总质量 kg
W_C=(w_C1)*M_TOTAL;                                     % 压缩过程压缩机总耗功 J

n_P=0.85;                                               % 水泵效率
syms s;                                                 % 压差控制
W_P=1/n_P*symsum((P_xun_in*Vs*((log((P_6+(s-1)*d)/P_xun_in))-(1-P_xun_in/(P_6+(s-1)*d))))+...
    V2*(P_6+s*d)*((log((P_6+s*d)/(P_6+(s-1)*d)))-(1-((P_6+(s-1)*d)/(P_6+s*d)))*(P_xun_in/(P_6+(s-1)*d))),1,n);
W_P=double(W_P);                                        % 压缩过程水泵总耗功 J

W_CHU=W_C+W_P;                                          % 储能过程总耗功 J (压缩机+循环压缩装置)
Q_COOL_1=M_TOTAL*(k*Rg/(k-1))*(T_C1_out-T_xun_in);      % 过程一放热量 J
Q_COOL=Q_COOL_1;                                        % 总放热量 J
%% 释能过程
kn=1.4;                         %循环膨胀过程指数
n_hy=0.85;                                              % 水轮机效率
Vs_shi=Vs;                                              % 循环释能装置单次释气量 m^3

P_XUN_SHI_OUT=1000000;

P_X=P_8_end;                                            % 储气罐进气前压力 Pa
Ww_1=0;                                                 % 迭代用，进气过程总输出功 J
Ww_2=0;                                                 % 迭代用，膨胀过程总输出功 J
nx=0;                                                   % 迭代用，排气次数
M_shi=0;                                                % 迭代用，排气总质量 kg
Q_HEAT_1=0;                                             % 迭代用，总加热量 J
while(P_X>P_6)
    P_XUN_IN=P_X;                                           % 循环释能装置进气结束时压力 Pa
    a=+inf;                                                 % 解方程参数
    b=+inf;                                                 % 解方程参数
    c=-inf;                                                 % 解方程参数
    while(abs(a+b-c)>0.001)
        V_SHI_0=P_XUN_SHI_OUT/P_XUN_IN*Vs_shi/((P_XUN_SHI_OUT/P_XUN_IN)^((kn-1)/kn));% 循环释能装置初始进气量 m^3
        a=P_XUN_IN*V_SHI_0/Rg/Th;                % 循环释能装置进气结束时气体质量 kg
        b=P_XUN_IN*V2/Rg/Tamb;                   % 储气罐进气后气体质量 kg
        c=P_X*V2/Rg/Tamb;                        % 储气罐进气前气体质量 kg
        instead=P_XUN_IN;                        % 循环释能装置进气结束时压力 Pa
        if a+b>c
            P_XUN_IN=P_XUN_IN-0.2;
        else
            P_XUN_IN=P_XUN_IN+0.2;
        end
    end
    Ww_n_2=-instead^(1/kn)*V_SHI_0/kn*(kn/(kn-1))*(P_XUN_SHI_OUT^((kn-1)/kn)-instead^((kn-1)/kn))-...
        instead^(1/kn)*V_SHI_0/kn*P_XUN_SHI_OUT*kn*(P_XUN_SHI_OUT^(-1/kn)-instead^(-1/kn));            % 单次膨胀过程总输出功 J
    Ww_n_1=-P_X*V2*log(instead/P_X)-P_X*V2*P_XUN_SHI_OUT*(1/instead)+P_X*V2*P_XUN_SHI_OUT*(1/P_X);% 单次进气过程总输出功 J
    P_X=instead;                                                                                  % 下一次进气前储气室压力 Pa
    Ww_2=Ww_2+Ww_n_2;                                                                             % 迭代用，膨胀过程总输出功 J
    Ww_1=Ww_1+(Th/Tamb)*Ww_n_1;                                                                             % 迭代用，进气过程总输出功 J
    T_XUN_shi_OUT=(P_XUN_SHI_OUT/instead)^((kn-1)/kn)*Th;                                           % 膨胀结束时气体温度 K
    nx=nx+1;                                                                                      % 迭代用，排气次数
    m_shi=instead*V_SHI_0/Rg/Th;                                                                  % 单次排气质量 kg
    M_shi=M_shi+m_shi;                                                                            % 迭代用，排气总质量 kg
    Q_n=m_shi*(k*Rg/(k-1))*(Th-T_XUN_shi_OUT);                                                    % 单次加热量 J
    Q_HEAT_1=Q_HEAT_1+Q_n;                                                                        % 迭代用，总加热量1 J
end
W_hy=n_hy*(Ww_2+Ww_1);                                                                            % 释能过程水轮机总输出功 J
% 透平
n_Tur=0.85;                                                                                       % 透平效率
epsilon_TUR=1/(P_XUN_SHI_OUT/P_amb);                                                              % 膨胀比
T_TUR1_OUT=Th*(1-n_Tur*(1-epsilon_TUR^((k-1)/k)));                                                % 透平1出口温度
w_TUR1=n_Tur*k*Rg*Th/(k-1)*(1-epsilon_TUR^((k-1)/k));                                             % 单位质量透平1做功 J
W_TUR1=w_TUR1*M_shi;                                                                              % 透平1做功 J 
W_TUR=(W_TUR1);                                                                                   % 透平机组总输出功 J
Q_HEAT_2=M_shi*(k*Rg/(k-1))*(Th-Tamb);
Q_HEAT=(Q_HEAT_1+Q_HEAT_2);                                                                         % 总加热量 J
%% 分析
n_1=W_hy/W_P;                         %  水轮机/水泵
n_2=(W_TUR+W_hy)/W_CHU;               % （透平+水轮机）/（压缩机+水泵）
n_3=(W_TUR+W_hy)/(Q_HEAT);            % （透平+水轮机）/（补热量）
n_4=(W_TUR+W_hy)/(W_CHU+Q_HEAT);      % （透平+水轮机）/（压缩机+水泵+补热量）
n_5=Q_COOL/W_CHU;                     %  放热量/（压缩机+水泵）
EXT_1=T_C1_out;                       %  压缩机1排气温度 K
EXT_3=T_TUR1_OUT;                     %  透平1排气温度 K

H00(i,j)=P_6;
H0(i,j)=P_8;
H1(i,j)=n_1;
H2(i,j)=n_2;
H3(i,j)=n_3;
H4(i,j)=n_4;
H5(i,j)=n_5;
H6(i,j)=EXT_1;
H8(i,j)=EXT_3;


% H0(i)=kn;
% H1(i)=n_1;
% H2(i)=n_2;
% H3(i)=n_3;
% H4(i)=n_4;
% H5(i)=n_5;
% H6(i)=EXT_1;
% H8(i)=EXT_3;
% HX=[H0;H1;H2;H3;H4;H5;H6;H8];
i,j
    end
end
%xlswrite('D:\DATE.xlsx',HX)
% fprintf('水轮机/水泵%2.2f\n',n_1)
% fprintf('（透平+水轮机）/（压缩机+水泵）%2.2f\n',n_2)
% fprintf('（透平+水轮机）/（补热量）%2.2f\n',n_3)
% fprintf('（透平+水轮机）/（压缩机+水泵+补热量）%2.2f\n',n_4)
% fprintf('放热量/（压缩机+水泵）%2.2f\n',n_5)
% fprintf('压缩机1排气温度为%2.2fK\n',EXT_1)
% fprintf('透平1排气温度为%2.2fK\n',EXT_3)
