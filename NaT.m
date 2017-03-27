function [gna] = NaT(Vm,dt,NumChan,stochastic)

t=dt:dt:(length(Vm).*dt);

% NA ARRAYS
m  = zeros(1,length(t));
h  = zeros(1,length(t));
m1h1 = zeros(1,length(t));
m2h1 = zeros(1,length(t));
m3h1 = zeros(1,length(t));
m4h1 = zeros(1,length(t));
m1h2 = zeros(1,length(t));
m2h2 = zeros(1,length(t));
m3h2 = zeros(1,length(t));
m4h2 = zeros(1,length(t));

gnasingle   = 0.000025; % ~20 pS is the estimated single channel conductance of NaT, NaP and DrK 
gnabar      = gnasingle.*NumChan;
N_NA        = NumChan;

i=1;

Vrest = Vm(1,i);
% Possible MATLAB equivalent expressions:
alpm(1,1) = 3.42.*( 0.11111.*(Vrest+33))./(1-exp(-0.11111.*(Vrest+33)));         
betm(1,1) = 27.6.*(-0.083333.*(Vrest+58))./(1-exp(0.083333.*(Vrest+58)));
beth(1,1) = 0.45.*( 0.11111.*(Vrest+21))./(1-exp(-0.11111.*(Vrest+21)));
alph(1,1) = 0.36.*(-0.083333.*(Vrest+48))./(1-exp(0.083333.*(Vrest+48)));

m(1,i)	= alpm(1,1) ./ (alpm(1,1) + betm(1,1));
h(1,i)	= alph(1,1) ./ (alph(1,1) + beth(1,1));

if stochastic==1
    m2 = round( [3.*alpm(1,1) ./ (3.*alpm(1,1) + betm(1,1))] .* N_NA );
    m3 = round( [2.*alpm(1,1) ./ (2.*alpm(1,1) + 2.*betm(1,1))] .* N_NA );
    m4 = round( [1.*alpm(1,1) ./ (1.*alpm(1,1) + 3.*betm(1,1))] .* N_NA );
    m1 = N_NA - [m2 + m3 + m4];
    m2h1(1,i) = round( m2 - m2 .* [alph(1,1) ./ (alph(1,1) + beth(1,1))] );
    m3h1(1,i) = round( m3 - m3 .* [alph(1,1) ./ (alph(1,1) + beth(1,1))] );
    m4h1(1,i) = round( m4 - m4 .* [alph(1,1) ./ (alph(1,1) + beth(1,1))] );
    m1h1(1,i) = round( m1 - m1 .* [alph(1,1) ./ (alph(1,1) + beth(1,1))] );
    m2h2(1,i) = m2 - m2h1(1,i);
    m3h2(1,i) = m3 - m3h1(1,i);
    m4h2(1,i) = m4 - m4h1(1,i);
    m1h2(1,i) = m1 - m1h1(1,i);
    N_Occup_Na = [m1h1(1,i),m2h1(1,i),m3h1(1,i),m4h1(1,i);m1h2(1,i),m2h2(1,i),m3h2(1,i),m4h2(1,i)];
end

for i = 2:length(Vm)

    
    % Possible MATLAB equivalent expressions:
    alpm(1,1) = 3.42.*(0.11111.*(Vm(1,i-1)+33))./(1-exp(-0.11111.*(Vm(1,i-1)+33)));         
    betm(1,1) = 27.6.*(-0.083333.*(Vm(1,i-1)+58))./(1-exp(0.083333.*(Vm(1,i-1)+58)));
    beth(1,1) = 0.45.*(0.11111.*(Vm(1,i-1)+21))./(1-exp(-0.11111.*(Vm(1,i-1)+21)));
    alph(1,1) = 0.36.*(-0.083333.*(Vm(1,i-1)+48))./(1-exp(0.083333.*(Vm(1,i-1)+48)));
        
    if stochastic == 1

        N_Occup_Na = [m1h1(i-1),m2h1(i-1),m3h1(i-1),m4h1(i-1);m1h2(i-1),m2h2(i-1),m3h2(i-1),m4h2(i-1)];

        P_Trans_m = [dt.*alpm(1,1).*[3,2,1];dt.*betm(1,1).*[1,2,3]];
        P_Trans_h = [dt.*alph(1,1);dt.*beth(1,1)];

        if max(N_Occup_Na(1,:))>N_NA || max(N_Occup_Na(2,:))>N_NA || min(N_Occup_Na(1,:))<0 || min(N_Occup_Na(2,:))<0
            disp('ERROR: Na Channel occupancies outside acceptable range');
            disp(N_Occup_Na);
            disp(i);
            return
        end

        for j=1:4

            if N_Occup_Na(1,j) > 0 
                tmp1 = rand(1,N_Occup_Na(1,j));
            else % Make sure no transitions can occur out of this state if it is unpopulated
                tmp1 = 1e67;              
            end

            if N_Occup_Na(2,j) > 0
                tmp2 = rand(1,N_Occup_Na(2,j));
            else % Make sure no transitions can occur out of this state if it is unpopulated
                tmp2 = 1e67;
            end


            switch j
                case 1,

                    tmp_1h = (tmp1 <= P_Trans_h(2,1));
                       T_Change_h(1,j) = sum(tmp_1h,2);
                    tmp3 = (tmp1 <= [P_Trans_m(1,j)+P_Trans_h(2,1)] & tmp1 > P_Trans_h(2,1));
                        T_Change_a(1,j) = sum(tmp3,2);

                    tmp_2h = (tmp2 < P_Trans_h(1,1));
                        T_Change_h(2,j) = sum(tmp_2h,2);
                    tmp4 = (tmp2 < [P_Trans_m(1,j)+P_Trans_h(1,1)] & tmp2 > P_Trans_h(1,1));
                        T_Change_a(2,j) = sum(tmp4,2);

                case {2,3},

                    tmp_1h = (tmp1 <= P_Trans_h(2,1));
                       T_Change_h(1,j) = sum(tmp_1h,2);
                    tmp3 = (tmp1 <= [P_Trans_m(1,j)+P_Trans_h(2,1)] & tmp1 > P_Trans_h(2,1));
                        T_Change_a(1,j) = sum(tmp3,2);
                    tmp5 = (tmp1 <= [P_Trans_m(2,j-1)+P_Trans_h(2,1)+P_Trans_m(1,j)] & tmp1 > [P_Trans_m(1,j)+P_Trans_h(2,1)]);
                        T_Change_b(1,j) = sum(tmp5,2);

                    tmp_2h = (tmp2 <= P_Trans_h(1,1));
                        T_Change_h(2,j) = sum(tmp_2h,2);
                    tmp4 = (tmp2 <= [P_Trans_m(1,j)+P_Trans_h(1,1)] & tmp2 > P_Trans_h(1,1));
                        T_Change_a(2,j) = sum(tmp4,2);
                    tmp6 = (tmp2 <= [P_Trans_m(2,j-1)+P_Trans_h(1,1)+P_Trans_m(1,j)] & tmp2 > [P_Trans_m(1,j)+P_Trans_h(1,1)]);
                        T_Change_b(2,j) = sum(tmp6,2);                                 

                case 4,

                    tmp_1h = (tmp1 <= P_Trans_h(2,1));
                       T_Change_h(1,j) = sum(tmp_1h,2);
                    tmp5 = (tmp1 <= [P_Trans_m(2,j-1)+P_Trans_h(2,1)] & tmp1 > P_Trans_h(2,1));
                        T_Change_b(1,j) = sum(tmp5,2);

                    tmp_2h = (tmp2 <= P_Trans_h(1,1));
                        T_Change_h(2,j) = sum(tmp_2h,2);
                    tmp6 = (tmp2 <= [P_Trans_m(2,j-1)+P_Trans_h(1,1)] & tmp2 > P_Trans_h(1,1));
                        T_Change_b(2,j) = sum(tmp6,2);  

            end

        end

        % Update states with number of transitions in - out
        m1h1(1,i)= m1h1(1,i-1) + T_Change_h(2,1) + T_Change_b(1,2) - T_Change_h(1,1) - T_Change_a(1,1);
        m1h2(1,i)= m1h2(1,i-1) + T_Change_h(1,1) + T_Change_b(2,2) - T_Change_h(2,1) - T_Change_a(2,1);

        m2h1(1,i)= m2h1(1,i-1) + T_Change_a(1,1) + T_Change_b(1,3) + T_Change_h(2,2) - T_Change_b(1,2) - T_Change_a(1,2) - T_Change_h(1,2);
        m2h2(1,i)= m2h2(1,i-1) + T_Change_a(2,1) + T_Change_b(2,3) + T_Change_h(1,2) - T_Change_b(2,2) - T_Change_a(2,2) - T_Change_h(2,2);

        m3h1(1,i)= m3h1(1,i-1) + T_Change_a(1,2) + T_Change_b(1,4) + T_Change_h(2,3) - T_Change_b(1,3) - T_Change_a(1,3) - T_Change_h(1,3);
        m3h2(1,i)= m3h2(1,i-1) + T_Change_a(2,2) + T_Change_b(2,4) + T_Change_h(1,3) - T_Change_b(2,3) - T_Change_a(2,3) - T_Change_h(2,3);

        m4h1(1,i)= m4h1(1,i-1) + T_Change_h(2,4) + T_Change_a(1,3) - T_Change_h(1,4) - T_Change_b(1,4);
        m4h2(1,i)= m4h2(1,i-1) + T_Change_h(1,4) + T_Change_a(2,3) - T_Change_h(2,4) - T_Change_b(2,4);

        gna(1,i) = gnasingle .* m4h1(1,i);

    else

        m(1,i)	 = m(1,i-1) + ( ( (alpm(1,1) .* (1 - m(1,i-1))) - (betm(1,1) .* m(1,i-1)) ) .* dt );
        h(1,i)	 = h(1,i-1) + ( ( (alph(1,1) .* (1 - h(1,i-1))) - (beth(1,1) .* h(1,i-1)) ) .* dt );

        gna(1,i) = gnabar .* m(1,i-1).^3 .* h(1,i-1);

    end
    
end

plot(t, gna);