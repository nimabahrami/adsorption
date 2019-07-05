clc
format short g

%% Variables

p_eva = 1401; 
p_cond = 4498;
T_cw_user = input('Please enter cold water temperature: (celsius) \n');
T_hw_user = input('Please enter hot water temperature: (celsius) \n');
t_cw = T_cw_user + 273;
t_hw = T_hw_user + 273;
m_sg = 6.75;
c_sg = 924;
c_water = 4180;
s1 = py.CoolProp.CoolProp.PropsSI('S','P',p_eva,'T',t_cw,'Water');
s2 = py.CoolProp.CoolProp.PropsSI('S','P',p_cond,'T',t_hw,'Water');
G1 = py.CoolProp.CoolProp.PropsSI('G','P',p_eva,'T',t_cw,'Water');
G2 = py.CoolProp.CoolProp.PropsSI('G','P',p_cond,'T',t_hw,'Water');
q_st = t_hw*(s2-s1)+G2-G1+100000;
k_0i = 2*10^(-12); % is assumed in article
r = 461.52;
n_cycle = 72;

%% Parameters
q_ads = q_st;   % is assumed in article
q_des = q_st;   % is assumed in article

%% Main
% user selection
fprintf('------------------choosing state----------------------- \n')
fprintf('please choose your mode: \n')
user = strcmp(input('AD + Cooling ? (Y/N) \n \n'),'Y'); % user selects status type exactly 'Y'

% Evaluation

if user == 1
    p1 = p_eva;
    p2 = p_cond;
    p3 = p_cond;
    p4 = p_eva;
    t1 = t_cw;
    t2 = abs(q_st/(r*(log(p1/p2)+(q_st/(r*t1)))));
    t3 = t_hw;
    t4 = q_st/(r*(log(p3/p4)+(q_st/(r*t3))));

    X1 = p1*k_0i*exp(q_st/(r*t1));
    X2 = p2*k_0i*exp(q_st/(r*t2));
    X3 = p3*k_0i*exp(q_st/(r*t3));
    X4 = p4*k_0i*exp(q_st/(r*t4));
    
    m_water_final = (X2 - X3)*m_sg;
    q_12 = (X1 * m_sg * c_water + m_sg * c_sg)*(t2 - t1);
    q_23 = (m_sg * c_sg + ((X2 + X3)/2) * m_sg * c_water)*(t3 - t2) + (X2 - X3) * m_sg * q_des;
    q_bed = q_12 + q_23;
    SEC = ((X1 * m_sg * c_water)*(t2-t1) + (m_sg * c_sg + ((X2+X3)/2)*m_sg*c_water)*(t3-t2)+(X2-X3)*m_sg*q_des)...
        /(X2*m_sg - X3*m_sg);
    q_34 = (X3*m_sg*c_water+m_sg*c_sg)*(t3-t4);
    q_41 = (m_sg*c_sg+((X4+X1)/2)*m_sg*c_water)*(t4-t1)+(X1-X4)*m_sg*q_ads;
    q_cooling = q_34 + q_41;
    cop = q_41/q_bed;
    sdwp = n_cycle * m_water_final;
    
    fprintf('Provided cooling is: \n')
    disp(q_cooling*n_cycle)
    fprintf('amount of water produced is: \n')
    disp(sdwp)
    fprintf('Specific Energy Consumption is: \n')
    disp(SEC*n_cycle)
    fprintf('COP is: \n')
    disp(cop)
else
    if user == 0
        user1 = strcmp(input('Just AD? (Y/N)'),'Y');
        if user1 == 1
           
            t1 = t_cw;
            t2 = t_hw;
            t3prime = t_cw;
            p1 = p_eva;
            p2 = p_eva;
            p3prime = p2*exp((q_st/(r*t2))-(q_st/(r*t3)));
            

            X1 = p1*k_0i*exp(q_st/(r*t1));
            X2 = p2*k_0i*exp(q_st/(r*t2));
            X3prime = X2;
            m_water = (X1 - X2)*m_sg;
            q_bed = q_12;
            q_12 = ( m_sg * c_sg + ((X1+X2)/2)*m_sg*c_water)*(t2 - t1) ...
                +(X1-X2)*m_sg*q_des;
            SEC = ((m_sg*c_sg+((X1+X2)/2)*m_sg*c_water)*(t2-t1)+(X1-X2)*m_sg*q_des) ...
                /((X1-X2)*m_sg);
        else
            p1 = p_eva;
            p2 = p_cond;
            p3 = p_cond;
            p4 = p_eva;
            t1 = t_cw;
            t2 = t_cw;
            t3 = t_hw;
            t4 = t_hw;

            X1 = p1*k_0i*exp(q_st/(r*t1));
            X2 = p2*k_0i*exp(q_st/(r*t2));
            X3 = p3*k_0i*exp(q_st/(r*t3));
            X4 = p4*k_0i*exp(q_st/(r*t4));
            m_water = (X1 - X3)*m_sg;
            q_23 = (m_sg*c_sg+((X2+X3)/2)*m_sg*c_water)*(t3-t2)+(X2-X3)*m_sg*q_des;
            q_bed = q_23;
            SEC = ((m_sg*c_sg+((X2+X3)/2)*m_sg*c_water)*(t3-t2)+(X2-X3)*m_sg*q_des) ...
                /((X1-X3)*m_sg);
        end
    end
end

%% ploting

fprintf('------------------Plotting effect ----------------------- \n')
user = strcmp(input('Do you want to plot the potable water mass \br in different temperature ranges \n for the first state? \n (Y/N) \n \n'),'Y');
if user == 1
    t1 = linspace(t_cw,50+273,10);
    t3 = ones(1,numel(t1))*t_hw;
    t2 = ones(1,numel(t1));
    t4 = ones(1,numel(t1));
    X2 = ones(1,numel(t1));
    X3 = ones(1,numel(t1));
    m_water = ones(1,numel(t1));
    cop = ones(1,numel(t1));
   
    for i=1:numel(t1)
        t2(i) = abs((-1*q_st)/(r*(log(p1/p2)+(q_st/(r*t1(i))))));
        t4(i) = q_st/(r*(log(p3/p4)+(q_st/(r*t3(i)))));
        X1(i) = p1*k_0i*exp(q_st/(r*t1(i)));
        X2(i) = p2*k_0i*exp(q_st/(r*t2(i)));
        X3(i) = p3*k_0i*exp(q_st/(r*t3(i)));
        X4(i) = p4*k_0i*exp(q_st/(r*t4(i)));
        m_water(i) = (X2(i) - X3(i))*m_sg*n_cycle;
        q_41(i) = (m_sg*c_sg+((X4(i)+X1(i))/2)*m_sg*c_water)*(t4(i)-t1(i))+(X1(i)-X4(i))*m_sg*q_ads;
        q_12(i) = (X1(i) * m_sg * c_water + m_sg * c_sg)*(t2(i) - t1(i));
        q_23(i) = (m_sg * c_sg + ((X2(i) + X3(i))/2) * m_sg * c_water)*(t3(i) - t2(i)) + (X2(i) - X3(i)) * m_sg * q_des;
        
        cop(i)= q_41(i)/(q_12(i)+q_23(i));
        SEC(i) = ((X1(i) * m_sg * c_water)*(t2(i)-t1(i)) + (m_sg * c_sg + ((X2(i)+X3(i))/2)*m_sg*c_water)*(t3(i)-t2(i))+(X2(i)-X3(i))*m_sg*q_des)...
        /(X2(i)*m_sg - X3(i)*m_sg);
    end
    
    figure(1)
    subplot(2,2,1)
    plot(t1,m_water,'k--')
	axis([t1(1)-5 t1(end)+5 0.05 1000])
    xlabel('Cool water temperature (K)')
    ylabel('Potable water output (kg)')
    title('the effect of increasing cold water temperature on mass of potable water')
    
    subplot(2,2,2)
    plot(t1,X2,'r--')
	axis([t1(1)-5 t1(end)+5 0 0.3])
    xlabel('Cool water temperature (K)')
    ylabel('amount adsorbed by adsbnt(kg/kg dry adsbnt)')
    title('the effect of increasing cold water temperature on X ')
    
    subplot(2,2,3)
    plot(t1,cop,'g--')
	axis([t1(1)-5 t1(end)+5 0.5 1])
    xlabel('Cool water temperature (K)')
    ylabel('COP')
    title('the effect of increasing cold water temperature on COP ')
    
    subplot(2,2,4)
    plot(t1,SEC,'g+')
	axis([t1(1)-5 t1(end)+5 2.7*10^6 4*10^6])
    xlabel('Cool water temperature (K)')
    ylabel('SEC')
    title('the effect of increasing cold water temperature on SEC ')
            
else
    fprintf('Okay! Good luck!')
end

%% comparison plots
% comparing the SEC and Cop and Mass for different cases
    
    m_water_1 = ones(1,10);
    t1_2 = linspace(t_cw,50+273,10);
    t2_2 = ones(1,numel(t1))*t_hw;
    t3_2 = ones(1,numel(t1))*t_cw;
    t4_2 = ones(1,numel(t1));
    p1_2 = ones(1,numel(t1))*p_eva;
    p2_2 = ones(1,numel(t1))*p_eva;
    p3_2 = ones(1,numel(t1))*p_eva;
    p4_2 = ones(1,numel(t1));
    X1_2 = ones(1,numel(t1));
    X2_2 = ones(1,numel(t1));
    X3_2 = ones(1,numel(t1));
    X4_2 = ones(1,numel(t1));
    m_water_2 = ones(1,numel(t1));
    SEC_2 = ones(1,numel(t1));
    
    t1_3 = linspace(t_cw,50+273,10);
    t3_3 = ones(1,numel(t1))*t_hw;
    t2_3 = ones(1,numel(t1))*t_cw;
    t4_3 = ones(1,numel(t1));
    p1_3 = ones(1,numel(t1))*p_eva;
    p2_3 = ones(1,numel(t1))*p_cond;
    p3_3 = ones(1,numel(t1))*p_cond;
    p4_3 = ones(1,numel(t1));
    X1_3 = ones(1,numel(t1));
    X2_3 = ones(1,numel(t1));
    X3_3 = ones(1,numel(t1));
    X4_3 = ones(1,numel(t1));
    m_water_3 = ones(1,numel(t1));
    SEC_3 = ones(1,numel(t1));
    
    for i=1:numel(t1)
        X1_2(i) = p1_2(i)*k_0i*exp(q_st/(r*t1_2(i)));
        X1_3(i) = p1_3(i)*k_0i*exp(q_st/(r*t1_3(i)));
        
        X2_2(i) = p2_2(i)*k_0i*exp(q_st/(r*t2_2(i)));
        X2_3(i) = p2_3(i)*k_0i*exp(q_st/(r*t2_3(i)));
        
        X3_2(i) = p3_2(i)*k_0i*exp(q_st/(r*t3_2(i)));
        X3_3(i) = p3_3(i)*k_0i*exp(q_st/(r*t3_3(i)))/10;
        
        m_water_1(i) = (X2(i) - X3(i))*m_sg;
        m_water_2(i) = (X1_2(i) - X2_2(i))*m_sg;
        m_water_3(i) = (X1_3(i) - X3_3(i))*m_sg +0.2 ;
    end

    figure(2)
    plot(t1,m_water_1,'r',t1,m_water_2,'g',t1,m_water_3,'b'),legend('1','2','3')
	axis([t1(1)-5 t1(end)+5 0.05 2.5])
    xlabel('Cool water temperature (K)')
    ylabel('Potable water output (kg)')
    title('comparison of output mass of potable water')

%% Cost evaluation

%data input
n = 30;                                  %plant lifetime
i_new = 0.05;                                %interest rate
j = 0;                                   %inflation
CRF = (i_new*(i_new+1)^n)/((i_new+1)^n-1);           %capital recovery factor
IWF = CRF;                               %inflation weighted factor= CRF because j=0

%amount of produced water
SDWP = m_water_final*72;                           %daily water production [m3/ton of adsorbent]
m_sg_1 = 6.75;                             %mass of silica gel per bed [kg]
m_sg_ton = (m_sg_1)/1000;                  %mass of siliga gel for 4 beds [ton]
Vl = SDWP*n*365*m_sg_ton;                    %amount of produced water for total lifetime of plant

%capital cost
Di1=43564516;                          %capital investment
Di2=Di1*CRF;
C_capital=Di2/Vl;                      % [$/m3]

%cost of electricity
Elec_cost=0.133;                       %cost of electricity
W_pump=1541.3;                         %kw
E_pump=W_pump*Elec_cost;               %cost of electricity per hour
C_electricity=(E_pump*8760)/Vl;        % [$/m3]

%cost of labor
a=0.0811;                              %percentage of direct capital investment
C_labor=(Di1*a*CRF)/Vl;                % [$/m3]

%maintenance cost
b=0.0463;                              %percentage of direct capital investment
C_maintenance=(Di1*b*CRF)/Vl;          % [$/m3]

%Results
LCOW=(C_capital+C_electricity+C_labor+C_maintenance);
fprintf('the levelized cost is: \n')
disp(LCOW)

%Error
Error=(LCOW-0.457)/0.457;


%% cost plot

SDWP= ones(1,numel(m_water));
V1 = ones(1,numel(m_water));
C_capital= ones(1,numel(m_water));
C_electricity= ones(1,numel(m_water));
C_labor= ones(1,numel(m_water));
C_maintenance= ones(1,numel(m_water));
LCOW= ones(1,numel(m_water));

for i=1:numel(m_water)
    SDWP(i) = m_water(i)*72;
    Vl(i) = SDWP(i)*n*365*m_sg;
    C_capital(i)=Di2/Vl(i);
    C_electricity(i)=(E_pump*8760)/Vl(i);
    C_labor(i)=(Di1*a*CRF)/Vl(i); 
    C_maintenance(i)=(Di1*b*CRF)/Vl(i);
    LCOW(i)=(C_capital(i)+C_electricity(i)+C_labor(i)+C_maintenance(i));
end

figure(3)
plot(m_water,LCOW)
xlabel('mass of water')
ylabel('LCOW')


