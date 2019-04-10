clear all
clc
format short g

%% Variables

p_eva = 1401;
p_cond = 4498;
t_cw = 25 + 273;
t_hw = 85 + 273;
m_sg = 6.75;
c_sg = 924;
c_water = 4180;
q_st = 2510000;
k_0i = 2*10^(-12);
r = 461.52;

%% Parameters
q_ads = q_st;
q_des = q_st;

%% Main
% user selection

fprintf('please choose your mode: \n')
user = strcmp(input('AD + Cooling ? (Y/N) \n \n'),'Y');

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
 publish('thermo.m','evalCode',false)        