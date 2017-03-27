function [currents,voltage] = GetSSIV(Vmin,Vmax,dV,conductanceVec);

% Define the voltage command
Vc = Vmin:dV:Vmax;

% Setting the cell's parameters for the simulation
Amemb   = 2.500005e-03 .* 2.500005e-03 .* pi .* 4;  % 50 um diameter
Cm		= Amemb.*1.67e03;                           % nF/cm^2

Vk		= -85;
Vna		=  55;
Vl		= -85;
Vhcn    = -30;
Vkca    = -85;

% Get the passed conductaces for the model
gNaT    = conductanceVec(1,1);
gNaP    = conductanceVec(1,2);
gKdr    = conductanceVec(1,3);
gKaf    = conductanceVec(1,4);
gKas    = conductanceVec(1,5);
gHcnWT  = conductanceVec(1,6);
gHcnKO  = conductanceVec(1,7);
gKca    = conductanceVec(1,8);
gLeak   = conductanceVec(1,9);

gnabar	= Amemb.*gNaT;         % about 6 times greater than Hoffman
gnapbar = Amemb.*gNaP;         % adding a persistent sodium current
gkbar	= Amemb.*gKdr;         % Same as Hoffman (all expressed in uS / cm^2)
gkafbar = Amemb.*gKaf;
gkasbar = Amemb.*gKas;
ghbar_wt= Amemb.*gHcnWT;
ghbar_ko= Amemb.*gHcnKO;
gkcabar = Amemb.*gKca;
glbar	= Amemb.*gLeak;

for i=1:length(Vc)

    % Possible MATLAB equivalent expressions:
    alpm(1,1) = 3.42.*( 0.11111.*(Vc(1,i)+33))./(1-exp(-0.11111.*(Vc(1,i)+33)));         
    betm(1,1) = 27.6.*(-0.083333.*(Vc(1,i)+58))./(1-exp(0.083333.*(Vc(1,i)+58)));
    beth(1,1) = 0.45.*( 0.11111.*(Vc(1,i)+21))./(1-exp(-0.11111.*(Vc(1,i)+21)));
    alph(1,1) = 0.36.*(-0.083333.*(Vc(1,i)+48))./(1-exp(0.083333.*(Vc(1,i)+48)));

    % Description from Matt's NEURON model using possible MATLAB equivalent expressions:        
    alpha_nap(1,1)  = 1.5./(1+exp(0.33333.*(-42.1-Vc(1,i))));
    beta_nap(1,1)   = 1./(1+exp(-0.3333.*(42.1-Vc(1,i))));
    beta_napi(1,1)  = 3.*0.000054757.*(0.38023.*(Vc(1,i)+64.409))./(1-exp(-0.38023.*(Vc(1,i)+64.409)));
    alpha_napi(1,1) = 3.*0.000040032.*(-0.21598.*(Vc(1,i)+17.014))./(1-exp(0.21598.*(Vc(1,i)+17.014)));

    % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
    alpn(1,1) = 0.2.*(0.1.*(Vc(1,i)+38))./(1-exp(-0.1.*(Vc(1,i)+38)));
    betn(1,1) = 0.6294.*(-0.02857.*(Vc(1,i)+47))./(1-exp(+0.02857.*(Vc(1,i)+47)));

    % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
    alpha_kaf(1,1) = 0.15.*(0.066667.*(Vc(1,i)+18.3))./(1-exp(-0.066667.*(Vc(1,i)+18.3)));
    beta_kaf(1,1) = 0.15.*(-0.066667.*(Vc(1,i)+18.3))./(1-exp(0.066667.*(Vc(1,i)+18.3)));
    alph_kaf(1,1) = 0.082.*(-0.121951.*(Vc(1,i)+58))./(1-exp(0.121951.*(Vc(1,i)+58)));
    beth_kaf(1,1) = 0.082.*(0.121951.*(Vc(1,i)+58))./(1-exp(-0.121951.*(Vc(1,i)+58)));            

   % The slowly inactivating potassium current
    alpha_kas(1,1) = alpha_kaf(1,1)./10;
    beta_kas(1,1)  = beta_kaf(1,1)./10;
    alph_kas(1,1)  = alph_kaf(1,1)./150;
    beth_kas(1,1)  = beth_kaf(1,1)./150;

    % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
    alpha_ko(1,1)   = 0.036./(1+exp(-0.044543.*(-148.7-Vc(1,i))));
    beta_ko(1,1)    = 0.0036./(1+exp(0.07728.*(-50.7-Vc(1,i))));           

    % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
    alpha_wt(1,1)   = 0.0366./(1+exp(-0.03.*(-118.75-Vc(1,i))));
    beta_wt(1,1)    = 0.066./(1+exp(0.06.*(-56.05-Vc(1,i))));
        
    % Description from Matt's NEURON model using possible MATLAB equivalent expressions:  
    alpha_kca(1,1)  = 0;
    beta_kca(1,1)   = 1.6;

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
        Inat(1,i)   = (gna(1,i)  .* (Vc(1,i) - Vna));
        Inap(1,i)   = (gnap(1,i) .* (Vc(1,i) - Vna));
        Idrk(1,i)   = (gk(1,i)   .* (Vc(1,i) - Vk));
        Ikaf(1,i)   = (gkaf(1,i) .* (Vc(1,i) - Vk));
        Ikas(1,i)   = (gkas(1,i) .* (Vc(1,i) - Vk));
        Ih(1,i)     = (ghcn(1,i) .* (Vc(1,i) - Vhcn));
        Ikca(1,i)   = (gkca(1,i) .* (Vc(1,i) - Vkca));
        Il(1,i)     = (gl(1,i)   .* (Vc(1,i) - Vl));

        Iion(1,i)	=  Inat(1,i) + Inap(1,i) + Idrk(1,i) + Ikaf(1,i) + Ikas(1,i) + Ih(1,i) + Ikca(1,i) + Il(1,i);
end
    
currents(1,:) = Inat(1,:);
currents(2,:) = Inap(1,:);
currents(3,:) = Idrk(1,:);
currents(4,:) = Ikas(1,:);
currents(5,:) = Ikaf(1,:);
currents(6,:) = Ih(1,:);
currents(7,:) = Ikca(1,:);
currents(8,:) = Il(1,:);
currents(9,:) = Iion(1,:);

voltage = Vc;
