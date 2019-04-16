clear all
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
q_st = 2510000; % is assumed in article
k_0i = 2*10^(-12); % is assumed in article
r = 461.52;

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
    
    m_water = (X2 - X3)*m_sg;
    q_12 = (X1 * m_sg * c_water + m_sg * c_sg)*(t2 - t1);
    q_23 = (m_sg * c_sg + ((X2 + X3)/2) * m_sg * c_water)*(t3 - t2) + (X2 - X3) * m_sg * q_des;
    q_bed = q_12 + q_23;
    SEC = ((X1 * m_sg * c_water)*(t2-t1) + (m_sg * c_sg + ((X2+X3)/2)*m_sg*c_water)*(t3-t2)+(X2-X3)*m_sg*q_des)...
        /(X2*m_sg - X3*m_sg);
    q_34 = (X3*m_sg*c_water+m_sg*c_sg)*(t3-t4);
    q_41 = (m_sg*c_sg+((X4+X1)/2)*m_sg*c_water)*(t4-t1)+(X1-X4)*m_sg*q_ads;
    q_cooling = q_34 + q_41;
    
    fprintf('Provided cooling is: \n')
    disp(q_cooling)
    fprintf('amount of water produced is: \n')
    disp(m_water)
    fprintf('Specific Energy Consumption is: \n')
    disp(SEC)
    
else if user == 0
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
   
    for i=1:numel(t1)
        t2(i) = abs((-1*q_st)/(r*(log(p1/p2)+(q_st/(r*t1(i))))));
        t4(i) = q_st/(r*(log(p3/p4)+(q_st/(r*t3(i)))));
        X2(i) = p2*k_0i*exp(q_st/(r*t2(i)));
        X3(i) = p3*k_0i*exp(q_st/(r*t3(i)));
        m_water(i) = (X2(i) - X3(i))*m_sg;
    end
    
    subplot(2,1,1)
    plot(t1,m_water,'k--')
	axis([t1(1)-5 t1(end)+5 0.05 2.5])
    xlabel('Cool water temperature (K)')
    ylabel('Potable water output (kg)')
    title('the effect of increasing cold water temperature on mass of potable water')
    
    subplot(2,1,2)
    plot(t1,X2,'r--')
	axis([t1(1)-5 t1(end)+5 0 0.3])
    xlabel('Cool water temperature (K)')
    ylabel('amount adsorbed by adsbnt(kg/kg dry adsbnt)')
    title('the effect of increasing cold water temperature on X ')
            
else fprintf('Okay! Good luck!')
end
