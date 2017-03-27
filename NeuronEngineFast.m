function [response,currents] = NeuronEngineFast(Vss,Iinj,dt,stochasticVec,conductanceVec,hdf,basename);
% TITLE: NeuronEngine
% 
% OVERVIEW:
%     A single compartment model neuron with stochastic ion channels
%         copyright: Joshua Dudman, 2005-8
%         feedback:  j.dudman@gmail.com
%         major update in 2009 to allow large pre-allocated random array
%         for fast execution of code
% 
% GENERAL:
%     This model uses a Markov-type model to solve for stochastic channels or can be toggled to give
%     deterministic channels using the HH formalism.
% 
% USAGE:
%     [response,currents] = FinalNeuronEngine_Ia(Vss,Iinj,dt,stochasticVec,conductanceVec);
%
% NOTES:
%     > Iinj is the injected current. The user must specify the length of the
%           time step and dt for Iinj. 
%     > A sampling (dt) of <= 0.02 ms is recommended for accurate
%     calculations of the rates
%     > Integration is either Runge-Kutta (RKFlag=1) or Back Euler (RKFlag=0)
% 
% PASSED FLAGS:
%     stochasticVec
%     toggles the channel conductance method
%         1 -> Markov model of stochastic channels
%         0 -> Deterministic HH channels
% 
%     conductanceVec
%     default model runs with the following conductances:
%         24000       NaT     % No change, 73598 channels, 25 pS
%         150         NaP     % No change, 471 channels, 25 pS
%         11000       Kdr     % No change, 34557 channels, 25 pS
%         100         Ka      % 100, 6545 channels, 6 ps
%         500         Kd      % 500, 436 channels, 18 ps
%         300         HCN1    % No change, 18064 channels, 1 ps
%         0           HCN2    % In ko version I use 70, 5498 channels, 1ps
%         750         KCa     % Assume small single channel like an SK-channel
%         174.075     Kl      % 139.26 density should be increased to 174.075 in new formalism

RKflag = 1;
tic
% Begin initializing some of the various matrices (greatly improves performance)
sweepsize= size(Iinj);
sweeps  = sweepsize(1,1);
respC   = cell(1,sweeps);
condC   = cell(1,sweeps);
gateC   = cell(1,sweeps);
e		= sweepsize(1,2).*dt;
t 		= dt:dt:e;
Vk		= -85;
Vna		=  55;
Vl		= -85;
Vhcn    = -30;
Vkca    = -85;

% CONDUCTANCES
gNaT    = conductanceVec(1,1);
gNaP    = conductanceVec(1,2);
gKdr    = conductanceVec(1,3);
gKaf    = conductanceVec(1,4);
gKas    = conductanceVec(1,5);
gHcnWT  = conductanceVec(1,6);
gHcnKO  = conductanceVec(1,7);
gKca    = conductanceVec(1,8);
gLeak   = conductanceVec(1,9);

% CELL PARAMETERS
Amemb   = 2.500005e-03 .* 2.500005e-03 .* pi .* 4;  % 50 um diameter
Cm		= Amemb.*1.67e03;                           % nF/cm^2

% GENERAL ARRAYS
tmp1    = zeros(1,80000);
tmp2    = zeros(1,80000);
tmp3    = zeros(1,80000);
tmp4    = zeros(1,80000);
tmp5    = zeros(1,80000);
tmp6    = zeros(1,80000);

Vm		= zeros(1,length(t));
Iion	= zeros(1,length(t));
dVdt	= zeros(1,length(t));
Inat    = zeros(1,length(t));
Inap    = zeros(1,length(t));
Idrk    = zeros(1,length(t));
Ikaf    = zeros(1,length(t));
Ikas    = zeros(1,length(t));
Ih      = zeros(1,length(t));
Ikca    = zeros(1,length(t));
Il      = zeros(1,length(t));
disp('Allocated general arrays')

% K ARRAYS
n       = zeros(1,length(t));
n1      = zeros(1,length(t));
n2      = zeros(1,length(t));
n3      = zeros(1,length(t));
n4      = zeros(1,length(t));
n5      = zeros(1,length(t));

n_kca   = zeros(1,length(t));
c_ca_1  = zeros(1,length(t));
c_ca_2  = zeros(1,length(t));
c_ca_3  = zeros(1,length(t));
o_ca    = zeros(1,length(t));

n_kas   = zeros(1,length(t));
h_kas   = zeros(1,length(t));
c_kas   = zeros(1,length(t));
o_kas   = zeros(1,length(t));
ci_kas  = zeros(1,length(t));
oi_kas  = zeros(1,length(t));

n_kaf   = zeros(1,length(t));
h_kaf   = zeros(1,length(t));
c_kaf   = zeros(1,length(t));
o_kaf   = zeros(1,length(t));
ci_kaf  = zeros(1,length(t));
oi_kaf  = zeros(1,length(t));

nl      = zeros(1,length(t));
cl      = zeros(1,length(t));
ol      = zeros(1,length(t));

gk      = zeros(1,length(t));
gkas    = zeros(1,length(t));
gkaf    = zeros(1,length(t));
gkca    = zeros(1,length(t));
gl      = zeros(1,length(t));
disp('Allocated K channel arrays');

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

m_nap  = zeros(1,length(t));
h_nap  = zeros(1,length(t));
c_nap  = zeros(1,length(t));
o_nap  = zeros(1,length(t));
ci_nap  = zeros(1,length(t));
oi_nap  = zeros(1,length(t));

gna  = zeros(1,length(t));
gnap = zeros(1,length(t));
disp('Allocated Na channel arrays');

% HCN ARRAYS
o_wt = zeros(1,length(t));
o_ko = zeros(1,length(t));
c_wt = zeros(1,length(t));
c_ko = zeros(1,length(t));
n_wt = zeros(1,length(t));
n_ko = zeros(1,length(t));

ghcn = zeros(1,length(t));
disp('Allocated HCN channel arrays');

%%% Begin sweeps of injected current loop
for inj_index=1:sweeps
    
    % Initialize the random number generator
    rand('twister',sum(100*clock));
    
    % Initialize membrane potential to steady-state
    Vrest = Vss;
    Vm(1,1)	= Vss;
    
    gnabar	= Amemb.*gNaT;         
    gnapbar = Amemb.*gNaP;         
    gkbar	= Amemb.*gKdr;         
    gkafbar = Amemb.*gKaf;
    gkasbar = Amemb.*gKas;
    ghbar_wt= Amemb.*gHcnWT;
    ghbar_ko= Amemb.*gHcnKO;
    gkcabar = Amemb.*gKca;
    glbar	= Amemb.*gLeak;        
    
    gnasingle   = 0.000025; % ~20 pS is the estimated single channel conductance of NaT, NaP and DrK 
    gnapsingle  = 0.000025; % ~20 pS is the estimated single channel conductance of NaT, NaP and DrK
    gksingle    = 0.000025; % ~20 pS is the estimated single channel conductance of NaT, NaP and DrK
    gkafsingle  = 0.000006; % From Chen and Johnston data
    gkassingle  = 0.000018; % From Chen and Johnston data
    ghcnsingle  = 0.000001; % Consistent with DiFrancesco estimate of very small (1 pS) unitary conductance
%     ghcnsingle  = 0.000025; % For testing effects of altered single channel conductance
    gkcasingle  = 0.000025; % Consistent with current being apamin-insensitive (larger conductance K channel)
    glsingle    = 0.000080; % Consistent with measurements of KCNK2
    
    N_NA        = round(gnabar ./ gnasingle);
    N_NaP       = round(gnapbar ./ gnapsingle);
    N_K         = round(gkbar ./ gksingle);		
    N_KAF       = round(gkafbar ./ gkafsingle);
    N_KAS       = round(gkasbar ./ gkassingle);		
    N_HCN_WT    = round(ghbar_wt ./ ghcnsingle);
    N_HCN_KO    = round(ghbar_ko ./ ghcnsingle);
    N_KCA       = round(gkcabar ./ gkcasingle);
    N_L         = round(glbar ./ glsingle);
    
%     kca_tau                     = 25;
    kca_tau                     = 50; % testing a large value measured experimentally
    spike_counter          = zeros(1,1);
    spike_counter(1,1)  = -inf;
    sn                              = 1;

% Display some online information about the progress of the simulation
	waitbarmsg = sprintf('Running sweep %g...',inj_index);
    waith = waitbar(0,waitbarmsg);
    drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Preallocate a huge array of random numbers that is at least as big as
%    the highest individual particle count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

AllChans = [N_NA,N_NaP,N_K,N_KAF,N_KAS,N_HCN_WT,N_HCN_KO,N_KCA,N_L];
LargestChan = max(AllChans)
RandomArray = rand(1,LargestChan.*20);
LargestChan = LargestChan.*19;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Begin the loop for numerical integration of the membrane current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for i=1:length(t)

	    if round(i./(length(t)./10))==i./(length(t)./10)
            waitbar(i./length(t));
            drawnow;
            disp(round(i./(length(t)./10)));
        end
        
%%%%%%%%%%%%%% Main numerical integration loop (begins after intialization) %%%%%%%%%%%%%%

		if i>1

        % Possible MATLAB equivalent expressions:
        alpm(1,1) = 3.42.*(0.11111.*(Vm(1,i-1)+33))./(1-exp(-0.11111.*(Vm(1,i-1)+33)));         
        betm(1,1) = 27.6.*(-0.083333.*(Vm(1,i-1)+58))./(1-exp(0.083333.*(Vm(1,i-1)+58)));
        beth(1,1) = 0.45.*(0.11111.*(Vm(1,i-1)+21))./(1-exp(-0.11111.*(Vm(1,i-1)+21)));
        alph(1,1) = 0.36.*(-0.083333.*(Vm(1,i-1)+48))./(1-exp(0.083333.*(Vm(1,i-1)+48)));
    
        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:        
        alpha_nap(1,1)  = 1.5./(1+exp(0.33333.*(-42.1-Vm(1,i-1))));
        beta_nap(1,1)   = 1./(1+exp(-0.3333.*(42.1-Vm(1,i-1))));
        beta_napi(1,1)  = 3.*0.000054757.*(0.38023.*(Vm(1,i-1)+64.409))./(1-exp(-0.38023.*(Vm(1,i-1)+64.409)));
        alpha_napi(1,1) = 3.*0.000040032.*(-0.21598.*(Vm(1,i-1)+17.014))./(1-exp(0.21598.*(Vm(1,i-1)+17.014)));
             
        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        alpn(1,1) = 0.2.*(0.1.*(Vm(1,i-1)+38))./(1-exp(-0.1.*(Vm(1,i-1)+38)));
        betn(1,1) = 0.6294.*(-0.02857.*(Vm(1,i-1)+47))./(1-exp(+0.02857.*(Vm(1,i-1)+47)));
                                             
        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        alpha_kaf(1,1)  = 0.15.*(0.066667.*(Vm(1,i-1)+18.3))./(1-exp(-0.066667.*(Vm(1,i-1)+18.3)));
        beta_kaf(1,1)   = 0.15.*(-0.066667.*(Vm(1,i-1)+18.3))./(1-exp(0.066667.*(Vm(1,i-1)+18.3)));
        alph_kaf(1,1)   = 0.082.*(-0.121951.*(Vm(1,i-1)+58))./(1-exp(0.121951.*(Vm(1,i-1)+58)));
        beth_kaf(1,1)   = 0.082.*(0.121951.*(Vm(1,i-1)+58))./(1-exp(-0.121951.*(Vm(1,i-1)+58)));            

        % The slowly inactivating potassium current
        alpha_kas(1,1) = alpha_kaf(1,1)./10;
        beta_kas(1,1)  = beta_kaf(1,1)./10;
        alph_kas(1,1)  = alph_kaf(1,1)./150;
        beth_kas(1,1)  = beth_kaf(1,1)./150;

        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        alpha_ko(1,1)   = 0.036./(1+exp(-0.044543.*(-148.7-Vm(1,i-1))));
        beta_ko(1,1)    = 0.0036./(1+exp(0.07728.*(-50.7-Vm(1,i-1))));           

        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
%         alpha_wt(1,1)   = 0.0366./(1+exp(-0.04918839154.*(-114.2-Vm(1,i-1))));
%         beta_wt(1,1)    = 0.066./(1+exp(0.09140767824.*(-51.5-Vm(1,i-1))));
        alpha_wt(1,1)   = 0.0366./(1+exp(-0.03.*(-118.75-Vm(1,i-1))));
        beta_wt(1,1)    = 0.066./(1+exp(0.06.*(-56.05-Vm(1,i-1))));
        % New Matt uses -118.75; 0.04 ;; -56.05; -0.09
        
        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        if Vm(1,i-1) > 0 & t(1,i-1) > spike_counter(1,length(spike_counter))+2
            spike_counter(1,sn) = t(1,i-1);
            alpha_kca(1,1)      = 1.5 .* sum(exp( (spike_counter(1,:)-t(1,i)) ./ kca_tau ), 2);
            sn                  = sn + 1;
        else
            if sn>1
%                 how_far         = find(spike_counter > (t(1,i)-(kca_tau.*8)));
                alpha_kca(1,1)  = 1.5 .* sum(exp( (spike_counter(1,:)-t(1,i)) ./ kca_tau ), 2);
            end
        end
%         beta_kca(1,1)   = 0.03;
        beta_kca(1,1)   = 1.6; % Changed to correspond to ~17.5 ms tau of AHP recovery

        % New leak conductance with kinetics:  
        alphaL(1,1) = 4;
        betaL(1,1)  = 1;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN THE CALCULATION OF GATING VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The sodium channel: 4 x 2 states (m+1 * h+1)
% m1h1, m2h1, m3h1, m4h1, m1h0, m2h0, m3h0, m4h0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if stochasticVec(1,1) == 1

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
                        %tmp1 = rand(1,N_Occup_Na(1,j));
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp1    = RandomArray(1,Lindex:Lindex+N_Occup_Na(1,j)-1);
                    else % Make sure no transitions can occur out of this state if it is unpopulated
                        tmp1 = 1e67;              
                    end

                    if N_Occup_Na(2,j) > 0
                        %tmp2 = rand(1,N_Occup_Na(2,j));
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp2    = RandomArray(1,Lindex:Lindex+N_Occup_Na(2,j)-1);
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The persistent sodium channel: 4 states
    % c_nap, o_nap, i_nap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if stochasticVec(1,2) == 1

                 N_NaP_Occup = [c_nap(1,i-1),o_nap(1,i-1);ci_nap(1,i-1),oi_nap(1,i-1)]; % update occupancies

                if (sum(N_NaP_Occup(1,:),2) + sum(N_NaP_Occup(2,:),2)) >N_NaP | min(N_NaP_Occup) < 0
                    disp('ERROR: exceeded total NaP channel count');
                    N_NaP_Occup
                    N_NaP
                    break
                end

                P_Trans_nap = [dt.*alpha_nap(1,1),dt.*beta_napi(1,1);dt.*beta_nap(1,1),dt.*alpha_napi(1,1)];

                for j=1:2

                    if N_NaP_Occup(1,j) > 0 % Transitions can only occur in the event that there are channels
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp1    = RandomArray(1,Lindex:Lindex+N_NaP_Occup(1,j)-1);
%                        tmp1 = rand(1,N_NaP_Occup(1,j));
                    else
                        tmp1 = 1e67;
                    end

                    if N_NaP_Occup(2,j) > 0 % Transitions can only occur in the event that there are channels
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp2    = RandomArray(1,Lindex:Lindex+N_NaP_Occup(2,j)-1);
                        %tmp2 = rand(1,N_NaP_Occup(2,j));
                    else
                        tmp2 = 1e67;
                    end

                    switch j

                        case 1, 
                            tmp3 = (tmp1 <= P_Trans_nap(1,1));
                                T_Change_NaP(1,j) = sum(tmp3,2);
                            tmp4 = (tmp1 <= [P_Trans_nap(1,1)+P_Trans_nap(1,2)] & tmp1 > P_Trans_nap(1,1));
                                T_Change_NaP_h(1,j) = sum(tmp4,2);

                            tmp5 = (tmp2 <= P_Trans_nap(1,1));
                                T_Change_NaP(2,j) = sum(tmp5,2);
                            tmp6 = (tmp2 <= [P_Trans_nap(1,1)+P_Trans_nap(2,2)] & tmp2 > P_Trans_nap(1,1));
                                T_Change_NaP_h(2,j) = sum(tmp6,2);                           

                        case 2,
                            tmp3 = (tmp1 <= P_Trans_nap(2,1));
                                T_Change_NaP(1,j) = sum(tmp3,2);
                            tmp4 = (tmp1 <= [P_Trans_nap(2,1)+P_Trans_nap(1,2)] & tmp1 > P_Trans_nap(2,1));
                                T_Change_NaP_h(1,j) = sum(tmp4,2);

                            tmp5 = (tmp2 <= P_Trans_nap(2,1));
                                T_Change_NaP(2,j) = sum(tmp5,2);
                            tmp6 = (tmp2 <= [P_Trans_nap(2,1)+P_Trans_nap(2,2)] & tmp2 > P_Trans_nap(2,1));
                                T_Change_NaP_h(2,j) = sum(tmp6,2);        
                    end

                end

                % Update states with number of transitions in - out
                c_nap(1,i)  = c_nap(1,i-1)  + T_Change_NaP(1,2) + T_Change_NaP_h(2,1) - T_Change_NaP(1,1) - T_Change_NaP_h(1,1);
                o_nap(1,i)  = o_nap(1,i-1)  + T_Change_NaP(1,1) + T_Change_NaP_h(2,2) - T_Change_NaP(1,2) - T_Change_NaP_h(1,2);
                ci_nap(1,i) = ci_nap(1,i-1) + T_Change_NaP(2,2) + T_Change_NaP_h(1,1) - T_Change_NaP(2,1) - T_Change_NaP_h(2,1);
                oi_nap(1,i) = oi_nap(1,i-1) + T_Change_NaP(2,1) + T_Change_NaP_h(1,2) - T_Change_NaP(2,2) - T_Change_NaP_h(2,2);

                % Solve the conductance value        
                gnap(1,i) = gnapsingle .* o_nap(1,i);

            else

                m_nap(1,i)  = m_nap(1,i-1) + (((alpha_nap(1,1)  .* (1 - m_nap(1,i-1))) - (beta_nap(1,1)  .* m_nap(1,i-1))) .* dt);
                h_nap(1,i)  = h_nap(1,i-1) + (((alpha_napi(1,1) .* (1 - h_nap(1,i-1))) - (beta_napi(1,1) .* h_nap(1,i-1))) .* dt);

                gnap(1,i)   = gnapbar .* m_nap(1,i-1) .* h_nap(1,i-1);

            end            
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The potassium channel: 5 states (n+1) 
% n1, n2, n3, n4, n5	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if stochasticVec(1,3) == 1
        N_Occup = [n1(i-1),n2(i-1),n3(i-1),n4(i-1),n5(i-1)];
        P_Trans = [dt.*alpn(1,1).*[4,3,2,1];dt.*betn(1,1).*[1,2,3,4]];

        if max(N_Occup)>N_K
            disp('Error: exceeded total drK channel count');
            break
        end

        for j=1:5

            if N_Occup(1,j) > 0
                Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                tmp1    = RandomArray(1,Lindex:Lindex+N_Occup(1,j)-1);
                %tmp1 = rand(1,N_Occup(1,j));
            else
                tmp1 = 1e67;
            end

            switch j
                case 1,
                    tmp2 = (tmp1 <= P_Trans(1,j));
                        T_Change(1,j) = sum(tmp2,2);
                case {2,3,4},
                    tmp2 = (tmp1 <= P_Trans(1,j));
                        T_Change(1,j) = sum(tmp2,2);
                    tmp3 = (tmp1 <= [P_Trans(2,j-1)+P_Trans(1,j)] & tmp1 > P_Trans(1,j));  
                        T_Change(2,j-1) = sum(tmp3,2);
                case 5,
                    tmp3 = (tmp1 <= P_Trans(2,j-1));
                        T_Change(2,j-1) = sum(tmp3,2);             
            end

        end

        n1(1,i)= n1(1,i-1) + T_Change(2,1) - T_Change(1,1);
        n2(1,i)= n2(1,i-1) + T_Change(1,1) + T_Change(2,2) - T_Change(1,2) - T_Change(2,1);
        n3(1,i)= n3(1,i-1) + T_Change(1,2) + T_Change(2,3) - T_Change(1,3) - T_Change(2,2);
        n4(1,i)= n4(1,i-1) + T_Change(1,3) + T_Change(2,4) - T_Change(1,4) - T_Change(2,3);
        n5(1,i)= n5(1,i-1) + T_Change(1,4) - T_Change(2,4);

        gk(i) = gksingle .* n5(i);
    else

        n(1,i)	= n(1,i-1) + (((alpn(1,1) .* (1 - n(1,i-1))) - (betn(1,1) .* n(1,i-1))) .* dt);
        gk(1,i) = gkbar .* n(1,i-1).^4;

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The kaf channel: 4 states
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if stochasticVec(1,4) == 1  
               
        N_Occup_Kaf = [c_kaf(1,i-1),o_kaf(1,i-1);ci_kaf(1,i-1),oi_kaf(1,i-1)];

        P_Trans_Kaf = [dt.*alpha_kaf(1,1),dt.*beth_kaf(1,1);dt.*beta_kaf(1,1),dt.*alph_kaf(1,1)];

        if max(N_Occup_Kaf(1,:))>N_KAF
            error('******ERROR******* exceeded total channel count in KS');
        end


        for j=1:2

            if N_Occup_Kaf(1,j) > 0 % Transitions can only occur in the event that there are channels
                Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                tmp1    = RandomArray(1,Lindex:Lindex+N_Occup_Kaf(1,j)-1);
                %tmp1 = rand(1,N_Occup_Kaf(1,j));
            else
                tmp1 = 1e67;
            end

            if N_Occup_Kaf(2,j) > 0 % Transitions can only occur in the event that there are channels
                Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                tmp2    = RandomArray(1,Lindex:Lindex+N_Occup_Kaf(2,j)-1);
                %tmp2 = rand(1,N_Occup_Kaf(2,j));
            else
                tmp2 = 1e67;
            end

            switch j

                case 1, 
                    tmp3 = (tmp1 <= P_Trans_Kaf(1,1));
                        T_Change_kaf(1,j) = sum(tmp3,2);
                    tmp4 = (tmp1 <= (P_Trans_Kaf(1,1)+P_Trans_Kaf(1,2)) & tmp1 > P_Trans_Kaf(1,1));
                        T_Change_kaf_h(1,j) = sum(tmp4,2);

                    tmp5 = (tmp2 <= P_Trans_Kaf(1,1));
                        T_Change_kaf(2,j) = sum(tmp5,2);
                    tmp6 = (tmp2 <= (P_Trans_Kaf(1,1)+P_Trans_Kaf(2,2)) & tmp2 > P_Trans_Kaf(1,1));
                        T_Change_kaf_h(2,j) = sum(tmp6,2);                           

                case 2,
                    tmp3 = (tmp1 <= P_Trans_Kaf(2,1));
                        T_Change_kaf(1,j) = sum(tmp3,2);
                    tmp4 = (tmp1 <= (P_Trans_Kaf(2,1)+P_Trans_Kaf(1,2)) & tmp1 > P_Trans_Kaf(2,1));
                        T_Change_kaf_h(1,j) = sum(tmp4,2);

                    tmp5 = (tmp2 <= P_Trans_Kaf(2,1));
                        T_Change_kaf(2,j) = sum(tmp5,2);
                    tmp6 = (tmp2 <= (P_Trans_Kaf(2,1)+P_Trans_Kaf(2,2)) & tmp2 > P_Trans_Kaf(2,1));
                        T_Change_kaf_h(2,j) = sum(tmp6,2);        
            end

        end

            % Update states with number of transitions in - out
            c_kaf(1,i)  = c_kaf(1,i-1)  + T_Change_kaf(1,2) + T_Change_kaf_h(2,1) - T_Change_kaf(1,1) - T_Change_kaf_h(1,1);
            o_kaf(1,i)  = o_kaf(1,i-1)  + T_Change_kaf(1,1) + T_Change_kaf_h(2,2) - T_Change_kaf(1,2) - T_Change_kaf_h(1,2);
            ci_kaf(1,i) = ci_kaf(1,i-1) + T_Change_kaf(2,2) + T_Change_kaf_h(1,1) - T_Change_kaf(2,1) - T_Change_kaf_h(2,1);
            oi_kaf(1,i) = oi_kaf(1,i-1) + T_Change_kaf(2,1) + T_Change_kaf_h(1,2) - T_Change_kaf(2,2) - T_Change_kaf_h(2,2);

            gkaf(1,i) = gkafsingle .* o_kaf(1,i);

    else
            n_kaf(1,i)	= n_kaf(1,i-1) + (((alpha_kaf(1,1) .* (1 - n_kaf(1,i-1))) - (beta_kaf(1,1) .* n_kaf(1,i-1))) .* dt);
            h_kaf(1,i)	= h_kaf(1,i-1) + (((alph_kaf(1,1) .* (1 - h_kaf(1,i-1))) - (beth_kaf(1,1) .* h_kaf(1,i-1))) .* dt);

            gkaf(1,i)    = gkafbar .* n_kaf(1,i-1) .* h_kaf(1,i-1);
    end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The kas channel: 4 states
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if stochasticVec(1,5) == 1
                
                N_Occup_Kas = [c_kas(1,i-1),o_kas(1,i-1);ci_kas(1,i-1),oi_kas(1,i-1)];

                P_Trans_Kas = [dt.*alpha_kas(1,1),dt.*beth_kas(1,1);dt.*beta_kas(1,1),dt.*alph_kas(1,1)];

                if max(N_Occup_Kas(1,:))>N_KAS
                    error('******ERROR******* exceeded total channel count in KS');
                end


                for j=1:2

                    if N_Occup_Kas(1,j) > 0 % Transitions can only occur in the event that there are channels
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp1    = RandomArray(1,Lindex:Lindex+N_Occup_Kas(1,j)-1);                        
                        %tmp1 = rand(1,N_Occup_Kas(1,j));
                    else
                        tmp1 = 1e67;
                    end

                    if N_Occup_Kas(2,j) > 0 % Transitions can only occur in the event that there are channels
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp2    = RandomArray(1,Lindex:Lindex+N_Occup_Kas(2,j)-1);
                        %tmp2 = rand(1,N_Occup_Kas(2,j));
                    else
                        tmp2 = 1e67;
                    end

                    switch j

                        case 1, 
                            tmp3 = (tmp1 <= P_Trans_Kas(1,1));
                                T_Change_kas(1,j) = sum(tmp3,2);
                            tmp4 = (tmp1 <= (P_Trans_Kas(1,1)+P_Trans_Kas(1,2)) & tmp1 > P_Trans_Kas(1,1));
                                T_Change_kas_h(1,j) = sum(tmp4,2);

                            tmp5 = (tmp2 <= P_Trans_Kas(1,1));
                                T_Change_kas(2,j) = sum(tmp5,2);
                            tmp6 = (tmp2 <= (P_Trans_Kas(1,1)+P_Trans_Kas(2,2)) & tmp2 > P_Trans_Kas(1,1));
                                T_Change_kas_h(2,j) = sum(tmp6,2);                           

                        case 2,
                            tmp3 = (tmp1 <= P_Trans_Kas(2,1));
                                T_Change_kas(1,j) = sum(tmp3,2);
                            tmp4 = (tmp1 <= (P_Trans_Kas(2,1)+P_Trans_Kas(1,2)) & tmp1 > P_Trans_Kas(2,1));
                                T_Change_kas_h(1,j) = sum(tmp4,2);

                            tmp5 = (tmp2 <= P_Trans_Kas(2,1));
                                T_Change_kas(2,j) = sum(tmp5,2);
                            tmp6 = (tmp2 <= (P_Trans_Kas(2,1)+P_Trans_Kas(2,2)) & tmp2 > P_Trans_Kas(2,1));
                                T_Change_kas_h(2,j) = sum(tmp6,2);        
                    end

                end

                % Update states with number of transitions in - out
                c_kas(1,i)  = c_kas(1,i-1)  + T_Change_kas(1,2) + T_Change_kas_h(2,1) - T_Change_kas(1,1) - T_Change_kas_h(1,1);
                o_kas(1,i)  = o_kas(1,i-1)  + T_Change_kas(1,1) + T_Change_kas_h(2,2) - T_Change_kas(1,2) - T_Change_kas_h(1,2);
                ci_kas(1,i) = ci_kas(1,i-1) + T_Change_kas(2,2) + T_Change_kas_h(1,1) - T_Change_kas(2,1) - T_Change_kas_h(2,1);
                oi_kas(1,i) = oi_kas(1,i-1) + T_Change_kas(2,1) + T_Change_kas_h(1,2) - T_Change_kas(2,2) - T_Change_kas_h(2,2);

                gkas(1,i) = gkassingle .* o_kas(1,i);
                
            else
                
                n_kas(1,i)	= n_kas(1,i-1) + (((alpha_kas(1,1) .* (1 - n_kas(1,i-1))) - (beta_kas(1,1) .* n_kas(1,i-1))) .* dt);
                h_kas(1,i)	= h_kas(1,i-1) + (((alph_kas(1,1)  .* (1 - h_kas(1,i-1))) - (beth_kas(1,1) .* h_kas(1,i-1))) .* dt);

                gkas(1,i)    = gkasbar .* n_kas(1,i-1) .* h_kas(1,i-1);
              
            end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The hcn channel: 2 x 2 states
% c0, c1, o0, o1  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
            if stochasticVec(1,6) == 1 | stochasticVec(1,7) == 1
                

                N_Occup_h = [c_ko(i-1),o_ko(i-1);c_wt(i-1),o_wt(i-1)];

                P_Trans_KO = [dt.*alpha_ko(1,1);dt.*beta_ko(1,1)];
                P_Trans_WT = [dt.*alpha_wt(1,1);dt.*beta_wt(1,1)];

                if max(N_Occup_h(1,:))>N_HCN_KO || max(N_Occup_h(2,:))>N_HCN_WT
                    disp('ERROR: exceeded total HCN channel count');
                    break
                end


                for j=1:2

                    if N_Occup_h(1,j) > 0
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp1    = RandomArray(1,Lindex:Lindex+N_Occup_h(1,j)-1);
                        %tmp1 = rand(1,N_Occup_h(1,j));
                    else % Make sure no transitions can occur out of this state
                        tmp1 = 1;              
                    end

                    if N_Occup_h(2,j) > 0
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp2    = RandomArray(1,Lindex:Lindex+N_Occup_h(2,j)-1);
                        %tmp2 = rand(1,N_Occup_h(2,j));
                    else % Make sure no transitions can occur out of this state
                        tmp2 = 1;
                    end

                    if j<2 % Calculate forward transitions in m
                        tmp3 = (tmp1 <= P_Trans_KO(1,1));
                        T_Change_Open(1,j) = sum(tmp3,2);
                        tmp4 = (tmp2 <= P_Trans_WT(1,1));
                        T_Change_Open(2,j) = sum(tmp4,2);
                    end

                    if j>1 % Calculate backward transitions in m
                        tmp5 = (tmp1 <= P_Trans_KO(2,1));
                        T_Change_Open(1,j) = sum(tmp5,2);
                        tmp6 = (tmp2 <= P_Trans_WT(2,1));
                        T_Change_Open(2,j) = sum(tmp6,2);
                    end

                end

            % Update states with number of transitions in - out
                % 1,1
                c_ko(1,i)= c_ko(1,i-1) + T_Change_Open(1,2) - T_Change_Open(1,1);
                % 2,1
                c_wt(1,i)= c_wt(1,i-1) + T_Change_Open(2,2) - T_Change_Open(2,1);
                % 1,2
                o_ko(1,i) = o_ko(1,i-1) + T_Change_Open(1,1) - T_Change_Open(1,2);
                % 2,2
                o_wt(1,i) = o_wt(1,i-1) + T_Change_Open(2,1) - T_Change_Open(2,2);

                ghcn(1,i) = ghcnsingle .* ( o_wt(1,i)  +  o_ko(1,i) );

            else

                n_wt(1,i) = [n_wt(1,i-1) + (((alpha_wt(1,1) .* (1 - n_wt(1,i-1))) - (beta_wt(1,1) .* n_wt(1,i-1))) .* dt)];
                n_ko(1,i) = [n_ko(1,i-1) + (((alpha_ko(1,1) .* (1 - n_ko(1,i-1))) - (beta_ko(1,1) .* n_ko(1,i-1))) .* dt)];

                ghcn(1,i)= [ghbar_wt .* n_wt(1,i-1)] + [ghbar_ko .* n_ko(1,i-1)];

            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The AHP potassium current
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if stochasticVec(1,8) == 1
                
                N_Occup_Ca = [c_ca_1(i-1),c_ca_2(i-1),c_ca_3(i-1),o_ca(i-1)];
                P_Trans_Ca = [dt.*alpha_kca(1,1).*[3,2,1];dt.*beta_kca(1,1).*[1,2,3]];

                if max(N_Occup_Ca)>N_KCA
                    disp('Error: exceeded total KCa channel count');
                    break
                end

                for j=1:4

                    if N_Occup_Ca(1,j) > 0
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp1    = RandomArray(1,Lindex:Lindex+N_Occup_Ca(1,j)-1);
                        %tmp1 = rand(1,N_Occup_Ca(1,j));
                    else
                        tmp1 = 1e67;
                    end

                    switch j
                        case 1,
                            tmp2 = (tmp1 <= P_Trans_Ca(1,j));
                                T_Change_Ca(1,j) = sum(tmp2,2);
                        case {2,3},
                            tmp2 = (tmp1 <= P_Trans_Ca(1,j));
                                T_Change_Ca(1,j) = sum(tmp2,2);
                            tmp3 = (tmp1 <= [P_Trans_Ca(2,j-1)+P_Trans_Ca(1,j)] & tmp1 > P_Trans_Ca(1,j));  
                                T_Change_Ca(2,j-1) = sum(tmp3,2);
                        case 4,
                            tmp3 = (tmp1 <= P_Trans_Ca(2,j-1));
                                T_Change_Ca(2,j-1) = sum(tmp3,2);             
                    end

                end

                c_ca_1(1,i)= c_ca_1(1,i-1) + T_Change_Ca(2,1)                    - T_Change_Ca(1,1);
                c_ca_2(1,i)= c_ca_2(1,i-1) + T_Change_Ca(1,1) + T_Change_Ca(2,2) - T_Change_Ca(1,2) - T_Change_Ca(2,1);
                c_ca_3(1,i)= c_ca_3(1,i-1) + T_Change_Ca(1,2) + T_Change_Ca(2,3) - T_Change_Ca(1,3) - T_Change_Ca(2,2);
                o_ca(1,i)  = o_ca(1,i-1)   + T_Change_Ca(1,3)                    - T_Change_Ca(2,3);

                gkca(i) = gkcasingle .* o_ca(i);

            else

                n_kca(1,i)	= n_kca(1,i-1) + (((alpha_kca(1,1) .* (1 - n_kca(1,i-1))) - (beta_kca(1,1) .* n_kca(1,i-1))) .* dt);
                gkca(1,i)   = gkcabar .* n_kca(1,i).^3;

            end
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The LEAK potassium current
%  Modeled as a KCNK2 like leak conductance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if stochasticVec(1,9) == 1
               
                N_Occup_l = [cl(i-1),ol(i-1)];
                P_Trans_Open = [dt.*alphaL(1,1);dt.*betaL(1,1)];

                if max(N_Occup_l(1,:))>N_L
                    disp('ERROR: exceeded total Leak channel count');
                    break
                end


                for j=1:2

                    if N_Occup_l(1,j) > 0
                        Lindex  = round((rand(1,1).*(LargestChan-1))+1);
                        tmp1    = RandomArray(1,Lindex:Lindex+N_Occup_l(1,j)-1);
                        %tmp1 = rand(1,N_Occup_l(1,j));
                    else % Make sure no transitions can occur out of this state
                        tmp1 = 1;              
                    end

                    if j<2 % Calculate forward transitions in m
                        tmp3 = (tmp1 <= P_Trans_Open(1,1));
                        T_Change_Open(1,j) = sum(tmp3,2);
                    end

                    if j>1 % Calculate backward transitions in m
                        tmp5 = (tmp1 <= P_Trans_Open(2,1));
                        T_Change_Open(1,j) = sum(tmp5,2);
                    end

                end

                % Update states with number of transitions in - out
                cl(1,i) = cl(1,i-1) + T_Change_Open(1,2) - T_Change_Open(1,1);
                ol(1,i) = ol(1,i-1) + T_Change_Open(1,1) - T_Change_Open(1,2);

                gl(1,i) = glsingle .* ol(1,i);

            else

                nl(1,i) = nl(1,i-1) + (((alphaL(1,1) .* (1 - nl(1,i-1))) - (betaL(1,1) .* nl(1,i-1))) .* dt);

                gl(1,i) = glbar .* nl(1,i-1);

            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END gating variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			% Determining the dV/dt using the HH net current (Im = 0) equation
            Inat(1,i)   = (gna(1,i)  .* (Vm(1,i-1) - Vna));
            Inap(1,i)   = (gnap(1,i) .* (Vm(1,i-1) - Vna));
            Idrk(1,i)   = (gk(1,i)   .* (Vm(1,i-1) - Vk));
            Ikas(1,i)   = (gkas(1,i) .* (Vm(1,i-1) - Vk));
            Ikaf(1,i)   = (gkaf(1,i) .* (Vm(1,i-1) - Vk));
            Ih(1,i)     = (ghcn(1,i) .* (Vm(1,i-1) - Vhcn));
            Ikca(1,i)   = (gkca(1,i) .* (Vm(1,i-1) - Vkca));
            Il(1,i)     = (gl(1,i)   .* (Vm(1,i-1) - Vl));
            
			Iion(1,i)	=  Inat(1,i) + Inap(1,i) + Idrk(1,i) + Ikas(1,i) + Ikaf(1,i) + Ih(1,i) + Ikca(1,i) + Il(1,i);
            
			% Toggle the calculation of dVdt to use either back Euler or 4th Order Runge-Kutta
			if RKflag == 1
				% Inidividual RK components
				F1 = (Iinj(inj_index,i-1) - Iion(i-1)) ./ Cm;
                F4 = ((Iinj(inj_index,i) - Iion(i)) ./ Cm);
				F2 = (F1+F4)./2;
				F3 = (F2+F4)./2;
				% Weighted sum of the Runge-Kutta components (4th order)
				dVdt(1,i) = ((F1 + (2.*F2) + (2.*F3) + F4) ./ 6);    
            else % Backward euler method (standard and quite stable)
				dVdt(1,i) = ( Iinj(inj_index,i-1)./Cm ) - ( Iion(1,i-1)./Cm );
			end
            
			Vm(1,i) 	= Vm(1,i-1) + (dVdt(1,i) .* dt);

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial run to get steady state values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        else
            
        % Possible MATLAB equivalent expressions:
        alpm(1,1) = 3.42.*( 0.11111.*(Vrest+33))./(1-exp(-0.11111.*(Vrest+33)));         
        betm(1,1) = 27.6.*(-0.083333.*(Vrest+58))./(1-exp(0.083333.*(Vrest+58)));
        beth(1,1) = 0.45.*( 0.11111.*(Vrest+21))./(1-exp(-0.11111.*(Vrest+21)));
        alph(1,1) = 0.36.*(-0.083333.*(Vrest+48))./(1-exp(0.083333.*(Vrest+48)));
            
        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:        
        alpha_nap(1,1)  = 1.5./(1+exp(0.33333.*(-42.1-Vrest)));
        beta_nap(1,1)   = 1./(1+exp(-0.3333.*(42.1-Vrest)));
        beta_napi(1,1)  = 3.*0.000054757.*(0.38023.*(Vrest+64.409))./(1-exp(-0.38023.*(Vrest+64.409)));
        alpha_napi(1,1) = 3.*0.000040032.*(-0.21598.*(Vrest+17.014))./(1-exp(0.21598.*(Vrest+17.014)));
            
        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        alpn(1,1) = 0.2.*(0.1.*(Vrest+38))./(1-exp(-0.1.*(Vrest+38)));
        betn(1,1) = 0.6294.*(-0.02857.*(Vrest+47))./(1-exp(+0.02857.*(Vrest+47)));
                                             
        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        alpha_kaf(1,1) = 0.15.*(0.066667.*(Vrest+18.3))./(1-exp(-0.066667.*(Vrest+18.3)));
        beta_kaf(1,1) = 0.15.*(-0.066667.*(Vrest+18.3))./(1-exp(0.066667.*(Vrest+18.3)));
        alph_kaf(1,1) = 0.082.*(-0.121951.*(Vrest+58))./(1-exp(0.121951.*(Vrest+58)));
        beth_kaf(1,1) = 0.082.*(0.121951.*(Vrest+58))./(1-exp(-0.121951.*(Vrest+58)));            

       % The slowly inactivating potassium current
        alpha_kas(1,1) = alpha_kaf(1,1)./10;
        beta_kas(1,1)  = beta_kaf(1,1)./10;
        alph_kas(1,1)  = alph_kaf(1,1)./150;
        beth_kas(1,1)  = beth_kaf(1,1)./150;

        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        alpha_ko(1,1)   = 0.036./(1+exp(-0.044543.*(-148.7-Vrest)));
        beta_ko(1,1)    = 0.0036./(1+exp(0.07728.*(-50.7-Vrest)));           

        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        alpha_wt(1,1)   = 0.0366./(1+exp(-0.03.*(-118.75-Vrest)));
        beta_wt(1,1)    = 0.066./(1+exp(0.06.*(-56.05-Vrest)));
                  
        
        % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
        alpha_kca(1,1)  = 0;        
        beta_kca(1,1)   = 0.8;
        
        alphaL(1,1) = 4;
        betaL(1,1) = 1;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        %%% Assuming current Vm has been constant for infinite time %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
            m(1,i)	= alpm(1,1) ./ (alpm(1,1) + betm(1,1));
            h(1,i)	= alph(1,1) ./ (alph(1,1) + beth(1,1));      
        
            m_nap(1,i)    = alpha_nap(1,1) ./ [alpha_nap(1,1)+beta_nap(1,1)];
            h_nap(1,i)    = alpha_napi(1,1)./ [alpha_napi(1,1)+beta_napi(1,1)];
            
            n(1,i)	= alpn(1,1) ./ (alpn(1,1) + betn(1,1));
            
            n_kaf(1,i) = alpha_kaf(1,1) ./ (beta_kaf(1,1) + alpha_kaf(1,1));
            h_kaf(1,i) = alph_kaf(1,1) ./ (beth_kaf(1,1) + alph_kaf(1,1));
           
            n_kas(1,i) = alpha_kas(1,1) ./ (beta_kas(1,1) + alpha_kas(1,1));
            h_kas(1,i) = alph_kas(1,1) ./ (beth_kas(1,1) + alph_kas(1,1));
                           
            n_wt(1,i) = alpha_wt(1,1) ./ (alpha_wt(1,1) + beta_wt(1,1));            
            n_ko(1,i) = alpha_ko(1,1) ./ (alpha_ko(1,1) + beta_ko(1,1));     
             
            n_kca(1,i)  = alpha_kca(1,1) ./ (beta_kca(1,1) + alpha_kca(1,1));

            nl(1,i)     = alphaL(1,1) ./ (alphaL(1,1) + betaL(1,1));
            
        % Converting gating back to conductance
            gna(1,i)    = gnabar .* m(1,i).^3 .* h(1,i);
            gnap(1,i)   = gnapbar .* m_nap(1,i) .* h_nap(1,i);
            gk(1,i)     = gkbar .* n(1,i).^4; 
            gkaf(1,i)   = gkafbar .* n_kaf(1,i) .* h_kaf(1,i);     
            gkas(1,i)   = gkasbar .* n_kas(1,i) .* h_kas(1,i);     
            ghcn(1,i)   = (ghbar_wt .* n_wt(1,i)) + (ghbar_ko .* n_ko(1,i));
            gkca(1,i)   = gkcabar .* n_kca(1,i).^3; 
            gl(1,i)     = glbar .* nl(1,i);
        
        % Determining the dV/dt using the HH net current (Im = 0) equation
            Inat(1,i)   = (gna(1,i)  .* (Vm(1,1) - Vna));
            Inap(1,i)   = (gnap(1,i) .* (Vm(1,1) - Vna));
            Idrk(1,i)   = (gk(1,i)   .* (Vm(1,1) - Vk));
            Ikaf(1,i)   = (gkaf(1,i) .* (Vm(1,1) - Vk));
            Ikas(1,i)   = (gkas(1,i) .* (Vm(1,1) - Vk));
            Ih(1,i)     = (ghcn(1,i) .* (Vm(1,1) - Vhcn));
            Ikca(1,i)   = (gkca(1,i) .* (Vm(1,1) - Vkca));
            Il(1,i)     = (gl(1,i)   .* (Vm(1,1) - Vl));
            
			Iion(1,i)	=  Inat(1,i) + Inap(1,i) + Idrk(1,i) + Ikaf(1,i) + Ikas(1,i) + Ih(1,i) + Ikca(1,i) + Il(1,i);
            
            if stochasticVec(1,1)==1
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
            
    
            if stochasticVec(1,2)==1
                o_nap(1,i) = round((h_nap(1,i)).*m_nap(1,i).*N_NaP);
                c_nap(1,i) = round((h_nap(1,i)).*(1-m_nap(1,i)).*N_NaP);
                oi_nap(1,i)= floor((1-h_nap(1,i)).*m_nap(1,i).*N_NaP);
                ci_nap(1,i)= floor((1-h_nap(1,i)).*(1-m_nap(1,i)).*N_NaP);
                N_NaP_Occup = [c_nap(1,i),o_nap(1,i);ci_nap(1,i),oi_nap(1,i)];
            end
            
            if stochasticVec(1,3)==1
                closed  = (betn(1,1) ./ (alpn(1,1) + betn(1,1))).^4;
                open    = (alpn(1,1) ./ (alpn(1,1) + betn(1,1))).^4;
                n5(1,i) = round(open .* N_K);
                n2(1,i) = round( ( (alpn(1,1).*4) ./ (alpn(1,1).*4 + betn(1,1).*4)) .* closed .* N_K );
                n3(1,i) = round( ( (alpn(1,1).*3) ./ (alpn(1,1).*4 + betn(1,1).*4)) .* closed .* N_K );
                n4(1,i) = round( ( (alpn(1,1).*2) ./ (alpn(1,1).*4 + betn(1,1).*4)) .* closed .* N_K );
                n1(1,i) = N_K - (n5(1,i) + n4(1,i) + n3(1,i) + n2(1,i));
                N_Occup = [n1(i),n2(i),n3(i),n4(i),n5(i)];
            end
            
            if stochasticVec(1,4)==1
                m_kaf_inf(1,i) =   alpha_kaf(1,1) ./ [alpha_kaf(1,1)+beta_kaf(1,1)];
                h_kaf_inf(1,i) =   alph_kaf(1,1) ./ [alph_kaf(1,1)+beth_kaf(1,1)];
                o_kaf(1,i) = round((1-h_kaf_inf(1,i)).*m_kaf_inf(1,i).*N_KAF);
                c_kaf(1,i) = round((1-h_kaf_inf(1,i)).*(1-m_kaf_inf(1,i)).*N_KAF);
                oi_kaf(1,i)= round(h_kaf_inf(1,i).*m_kaf_inf(1,i).*N_KAF);
                ci_kaf(1,i)= round(h_kaf_inf(1,i).*(1-m_kaf_inf(1,i)).*N_KAF);
                N_Occup_Kaf = [c_kaf(1,i),o_kaf(1,i);ci_kaf(1,i),oi_kaf(1,i)];
            end
            
            if stochasticVec(1,5)==1
                m_kas_inf(1,i) =   alpha_kas(1,1) ./ [alpha_kas(1,1)+beta_kas(1,1)];
                h_kas_inf(1,i) =   alph_kas(1,1) ./ [alph_kas(1,1)+beth_kas(1,1)];
                o_kas(1,i) = round((h_kas_inf(1,i)).*m_kas_inf(1,i).*N_KAS);
                c_kas(1,i) = round((h_kas_inf(1,i)).*(1-m_kas_inf(1,i)).*N_KAS);
                oi_kas(1,i)= round((1-h_kas_inf(1,i)).*m_kas_inf(1,i).*N_KAS);
                ci_kas(1,i)= round((1-h_kas_inf(1,i)).*(1-m_kas_inf(1,i)).*N_KAS);
                N_Occup_Kas = [c_kas(1,i),o_kas(1,i);ci_kas(1,i),oi_kas(1,i)];
            end
            
            if stochasticVec(1,6)==1 | stochasticVec(1,7)==1
                o_wt(1,i) = round( n_wt(1,1).*N_HCN_WT );
                o_ko(1,i) = round( n_ko(1,1).*N_HCN_KO );
                c_wt(1,i) = round( N_HCN_WT - o_wt(1,i) );
                c_ko(1,i) = round( N_HCN_KO - o_ko(1,i) );    
                N_Occup_h = [c_ko(1,i),o_ko(1,i);c_wt(1,i),o_wt(1,i)];
            end
            
            if stochasticVec(1,8)==1
                c_ca_2(1,i)= round(N_KCA.*(n_kca(1,i).^1));
                c_ca_3(1,i)= round(N_KCA.*(n_kca(1,i).^2));
                o_ca(1,i)  = round(N_KCA.*(n_kca(1,i).^3));
                c_ca_1(1,i)= round(N_KCA - c_ca_2(1,i)+ c_ca_3(1,i) + o_ca(1,i));
                N_Occup_Ca = [c_ca_1(1,i),c_ca_2(1,i),c_ca_3(1,i),o_ca(1,i)];
            end
            
            if stochasticVec(1,9)==1
                cl(1,i) = round(0.2.*N_L);
                ol(1,i) = round(0.8.*N_L);
                N_Occup_l = [cl(1,i),ol(1,i)];
            end
            
		end
        
	end
    
    
    response(1,:) = Vm(1,:);
    response(2,:) = t(1,:);
    response(3,:) = Iinj(inj_index,:);

    currents(1,:) = Inat(1,:);
    currents(2,:) = Inap(1,:);
    currents(3,:) = Idrk(1,:);
    currents(4,:) = Ikas(1,:);
    currents(5,:) = Ikaf(1,:);
    currents(6,:) = Ih(1,:);
    currents(7,:) = Ikca(1,:);
    currents(8,:) = Il(1,:);
    currents(9,:) = Iion(1,:);

%     currents(1,1) = 0; % if required for memory space in long simulations

    disp(sprintf('Completed sweep %g of %g',inj_index,sweeps));
    
    disp(['Spike Count: ' num2str(sn-1)]);
%     disp(['Spike Times: ' num2str(spike_counter)]);

    if hdf
        if inj_index==1
            filenamev = sprintf('%s_v.h5',basename);
            filenamei = sprintf('%s_i.h5',basename);
            name = sprintf('/RecordA%g',inj_index);
            hdf5write(filenamev,name,response(1,:));
            name = sprintf('/RecordB%g',inj_index);
            hdf5write(filenamei,name,response(3,:));

        else
            filenamev = sprintf('%s_v.h5',basename);
            filenamei = sprintf('%s_i.h5',basename);
            name = sprintf('/RecordA%g',inj_index);
            hdf5write(filenamev,name,response(1,:),'WriteMode','append');
            name = sprintf('/RecordB%g',inj_index);
            hdf5write(filenamei,name,response(3,:),'WriteMode','append');

        end
    end

    eval(sprintf('save %s_%g response currents',basename,inj_index));

    close(waith);
            
    figure(1);
    plot(response(2,:),response(1,:));
    drawnow; hold on;
    
end
toc