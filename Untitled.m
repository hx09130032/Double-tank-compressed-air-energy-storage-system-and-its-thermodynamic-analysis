%2020.12.19
%% ������ѹ����+����ѭ��ѹ����/���������ͻ�+���ѭ�����ͣ�
clc;
clear all;
close all;
for i=1:11
    for j=1:11
%% �Է����ģ���������ݻ���С
        % ������ѹ�� 4-6MPa
        % �� 4MPa ѹ���� Th �¶Ƚ���1MW��ģ����1h(3.6*10^9 J)
        
k=1.4;                                       % ����ָ��
Rg=287;                                      % ���峣��
Tamb=273.15+20;                              % �����¶� K
P_amb=101325;                                % ����ѹ�� Pa
Th=273.15+300;                               % ̫���ܼ��Ⱥ������¶� K
n_T=0.85;                                    % ͸ƽЧ��
P_8=5500000+100000*(j-1);                                 % ���������ѹ�� Pa
P_6=3500000+100000*(i-1);                                 % ���������ѹ�� Pa
sss=sqrt(P_6/P_amb);
T_OUT=Th*(1-n_T*(1-(1/sss)^((k-1)/k)));     % ͸ƽ�����¶� K
w=k*Rg*Th/(k-1)*n_T*(1-(1/sss)^((k-1)/k));  % ��λ���������������� J
M=(3.6*10^9)/(2*w);                          % ���������� kg
m=M/3600;                                    % ������������ kg/s
         % P_8*V-P_6*V=M*Rg*Tamb
V=M*Rg*Tamb/(P_8-P_6);                       % �������ݻ� m^3
%% ���ܹ���
Vs=5;                                                   % ѭ��ѹ��װ�õ���ѹ���� m^3
V2=V;                                                   % �������ݻ� m^3
n_C=0.85;                                               % ѹ����Ч��
          % һ��ѹ��         
P_xun_in=400000;        

epsilon=(P_xun_in/P_amb);                               % ѹ��
T_C1_out=Tamb+[(epsilon)^((k-1)/k)-1]*Tamb/n_C;         % ��һ��ѹ���������¶� K

T_xun_in=Tamb;                                          % ѭ��ѹ��װ�ý����¶� K

d=P_xun_in*Tamb*Vs/T_xun_in/V2;                         % ѭ��ѹ��������ѹ�� Pa
P_8_end=(fix((P_8-P_6)/d)+1)*d+P_6;                     % ��������ʵ����ѹ�� Pa
n=(P_8_end-P_6)/d;                                      % ѹ������
w_C1=(1/n_C)*Rg*Tamb*(k/(k-1))*[epsilon^((k-1)/k)-1];   % ��λ����һ��ѹ���Ĺ� J
M_TOTAL=P_8_end*V2/Rg/Tamb-P_6*V2/Rg/Tamb;              % ѹ���������������� kg
W_C=(w_C1)*M_TOTAL;                                     % ѹ������ѹ�����ܺĹ� J

n_P=0.85;                                               % ˮ��Ч��
syms s;                                                 % ѹ�����
W_P=1/n_P*symsum((P_xun_in*Vs*((log((P_6+(s-1)*d)/P_xun_in))-(1-P_xun_in/(P_6+(s-1)*d))))+...
    V2*(P_6+s*d)*((log((P_6+s*d)/(P_6+(s-1)*d)))-(1-((P_6+(s-1)*d)/(P_6+s*d)))*(P_xun_in/(P_6+(s-1)*d))),1,n);
W_P=double(W_P);                                        % ѹ������ˮ���ܺĹ� J

W_CHU=W_C+W_P;                                          % ���ܹ����ܺĹ� J (ѹ����+ѭ��ѹ��װ��)
Q_COOL_1=M_TOTAL*(k*Rg/(k-1))*(T_C1_out-T_xun_in);      % ����һ������ J
Q_COOL=Q_COOL_1;                                        % �ܷ����� J
%% ���ܹ���
kn=1.4;                         %ѭ�����͹���ָ��
n_hy=0.85;                                              % ˮ�ֻ�Ч��
Vs_shi=Vs;                                              % ѭ������װ�õ��������� m^3

P_XUN_SHI_OUT=1000000;

P_X=P_8_end;                                            % �����޽���ǰѹ�� Pa
Ww_1=0;                                                 % �����ã���������������� J
Ww_2=0;                                                 % �����ã����͹���������� J
nx=0;                                                   % �����ã���������
M_shi=0;                                                % �����ã����������� kg
Q_HEAT_1=0;                                             % �����ã��ܼ����� J
while(P_X>P_6)
    P_XUN_IN=P_X;                                           % ѭ������װ�ý�������ʱѹ�� Pa
    a=+inf;                                                 % �ⷽ�̲���
    b=+inf;                                                 % �ⷽ�̲���
    c=-inf;                                                 % �ⷽ�̲���
    while(abs(a+b-c)>0.001)
        V_SHI_0=P_XUN_SHI_OUT/P_XUN_IN*Vs_shi/((P_XUN_SHI_OUT/P_XUN_IN)^((kn-1)/kn));% ѭ������װ�ó�ʼ������ m^3
        a=P_XUN_IN*V_SHI_0/Rg/Th;                % ѭ������װ�ý�������ʱ�������� kg
        b=P_XUN_IN*V2/Rg/Tamb;                   % �����޽������������� kg
        c=P_X*V2/Rg/Tamb;                        % �����޽���ǰ�������� kg
        instead=P_XUN_IN;                        % ѭ������װ�ý�������ʱѹ�� Pa
        if a+b>c
            P_XUN_IN=P_XUN_IN-0.2;
        else
            P_XUN_IN=P_XUN_IN+0.2;
        end
    end
    Ww_n_2=-instead^(1/kn)*V_SHI_0/kn*(kn/(kn-1))*(P_XUN_SHI_OUT^((kn-1)/kn)-instead^((kn-1)/kn))-...
        instead^(1/kn)*V_SHI_0/kn*P_XUN_SHI_OUT*kn*(P_XUN_SHI_OUT^(-1/kn)-instead^(-1/kn));            % �������͹���������� J
    Ww_n_1=-P_X*V2*log(instead/P_X)-P_X*V2*P_XUN_SHI_OUT*(1/instead)+P_X*V2*P_XUN_SHI_OUT*(1/P_X);% ���ν�������������� J
    P_X=instead;                                                                                  % ��һ�ν���ǰ������ѹ�� Pa
    Ww_2=Ww_2+Ww_n_2;                                                                             % �����ã����͹���������� J
    Ww_1=Ww_1+(Th/Tamb)*Ww_n_1;                                                                             % �����ã���������������� J
    T_XUN_shi_OUT=(P_XUN_SHI_OUT/instead)^((kn-1)/kn)*Th;                                           % ���ͽ���ʱ�����¶� K
    nx=nx+1;                                                                                      % �����ã���������
    m_shi=instead*V_SHI_0/Rg/Th;                                                                  % ������������ kg
    M_shi=M_shi+m_shi;                                                                            % �����ã����������� kg
    Q_n=m_shi*(k*Rg/(k-1))*(Th-T_XUN_shi_OUT);                                                    % ���μ����� J
    Q_HEAT_1=Q_HEAT_1+Q_n;                                                                        % �����ã��ܼ�����1 J
end
W_hy=n_hy*(Ww_2+Ww_1);                                                                            % ���ܹ���ˮ�ֻ�������� J
% ͸ƽ
n_Tur=0.85;                                                                                       % ͸ƽЧ��
epsilon_TUR=1/(P_XUN_SHI_OUT/P_amb);                                                              % ���ͱ�
T_TUR1_OUT=Th*(1-n_Tur*(1-epsilon_TUR^((k-1)/k)));                                                % ͸ƽ1�����¶�
w_TUR1=n_Tur*k*Rg*Th/(k-1)*(1-epsilon_TUR^((k-1)/k));                                             % ��λ����͸ƽ1���� J
W_TUR1=w_TUR1*M_shi;                                                                              % ͸ƽ1���� J 
W_TUR=(W_TUR1);                                                                                   % ͸ƽ����������� J
Q_HEAT_2=M_shi*(k*Rg/(k-1))*(Th-Tamb);
Q_HEAT=(Q_HEAT_1+Q_HEAT_2);                                                                         % �ܼ����� J
%% ����
n_1=W_hy/W_P;                         %  ˮ�ֻ�/ˮ��
n_2=(W_TUR+W_hy)/W_CHU;               % ��͸ƽ+ˮ�ֻ���/��ѹ����+ˮ�ã�
n_3=(W_TUR+W_hy)/(Q_HEAT);            % ��͸ƽ+ˮ�ֻ���/����������
n_4=(W_TUR+W_hy)/(W_CHU+Q_HEAT);      % ��͸ƽ+ˮ�ֻ���/��ѹ����+ˮ��+��������
n_5=Q_COOL/W_CHU;                     %  ������/��ѹ����+ˮ�ã�
EXT_1=T_C1_out;                       %  ѹ����1�����¶� K
EXT_3=T_TUR1_OUT;                     %  ͸ƽ1�����¶� K

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
% fprintf('ˮ�ֻ�/ˮ��%2.2f\n',n_1)
% fprintf('��͸ƽ+ˮ�ֻ���/��ѹ����+ˮ�ã�%2.2f\n',n_2)
% fprintf('��͸ƽ+ˮ�ֻ���/����������%2.2f\n',n_3)
% fprintf('��͸ƽ+ˮ�ֻ���/��ѹ����+ˮ��+��������%2.2f\n',n_4)
% fprintf('������/��ѹ����+ˮ�ã�%2.2f\n',n_5)
% fprintf('ѹ����1�����¶�Ϊ%2.2fK\n',EXT_1)
% fprintf('͸ƽ1�����¶�Ϊ%2.2fK\n',EXT_3)
