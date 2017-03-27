function [response,currents,Iinj] = RunControl(ss_flag,numsweeps,duration,on,off,amplitude,IinjType,stochasticVec,conductanceVec,hdf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% DETERMINE SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ss_flag

    [Vss] = LookupSSvoltage(conductanceVec,amplitude(1,1));
%     [Vss] = LookupSSvoltage(conductanceVec,amplitude(1,2));
    
else
    
    [Vss] = LookupSSvoltage(conductanceVec,0);
    
end
 
% Vss = -69.677 FROM MATT'S SIM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% TEST PARAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Default values for the conductance vector are ([NaT,NaP,Kdr,Ka,Kd,HCN1,HCN2,KCa,Kl]):
    % conductanceVec = [24000,150,11000,100,500,230, 0,100,174.075];
    % conductanceVec = [24000,150,11000,100,500,  0,70,100,174.075];
    % After fixing the AHP halfwidth:
    % conductanceVec = [24000,150,11000,100,500,250, 0,425,174.075];
    % conductanceVec = [24000,150,11000,100,500,  0,75,425,174.075];
    % Corrected conductances to match slope conductance data
    % conductanceVec = [24000,75,11000,100,500,600, 0,425,74.075];
    % conductanceVec = [24000,75,11000,100,500,0, 50,425,74.075];
    % stochasticVec  = [1,1,1,1,1,1,1,1,1];
    % stochasticVec  = [0,0,0,0,0,0,0,0,0];

%     conductanceVec = [24000,75,11000,100,500,  500, 0,425,150];
%     conductanceVec = [24000,75,11000,100,500,  0, 40,425,150];


    test2 = size(conductanceVec);
    if test2(1,2) ~= 9
        error('The conductance vector is the wrong size, please refer to CurrentDensities m-file')
    end

    test2 = size(stochasticVec);
    if test2(1,2) ~= 9
        error('The stochastic flag vector is the wrong size, please refer to CurrentDensities m-file')
    end

    if stochasticVec(1,1)==1 | stochasticVec(1,3)==1
        disp('Buyer beware: there are MANY NaT and Kdr channels in the model and the stochastic version runs quite slowly.')
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% DEFAULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt      = 0.01; % ms
    RKflag  = 1;    % use Runge-Kutta fourth order integration method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% CREATE STIMULUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % The amplitude argument can be a matrix in which the first column are DC
    % amplitudes and the second column is signal amplitude. The function can
    % generally tell the difference between 1,1 : 1,2 : num,2 matrices. In the
    % event of a num,1 or num,2 matrix delta is ignored for 'noisy'.

    test = size(amplitude);

    Iinj = zeros(numsweeps,duration./dt);

    switch lower(IinjType)
        
      case 'user'
          Iinj = amplitude;
          
      case 'zap'
        t = (on./dt):(off./dt);
        Itmp = amplitude(1,1) .* sin( (2 .* pi.* (t.*dt) .* (0.000004 .* (t.*dt))) + 0.0001);

        for i=1:numsweeps
%             Iinj(i,(on./dt):(off./dt)) = Itmp + (amplitude(1,2)*(i-1));
             Iinj(i,(on./dt):(off./dt)) = Itmp + (amplitude(1,2));
        end

      case 'dc'
        delta = amplitude(1,2);
        for i=1:numsweeps
            Iinj(i,(on./dt):(off./dt)) = amplitude(1,1) + (delta.*(i-1));
        end

      case 'ramp'
        begini  = amplitude(1,1);
        endi    = amplitude(1,2);
        slopei  = (endi-begini) ./ ((off./dt)-(on./dt));
        for i=1:numsweeps
            Iinj(i,(on./dt):(off./dt)) = slopei.*[0:((off./dt)-(on./dt))] + begini;
        end


      case 'noisy'
        for i=1:numsweeps
          if test(1,1) == numsweeps

            if test(1,2)==2
                Itmp = amplitude(i,2) .* randn(1,length(Iinj));
            else
                Itmp = amplitude(i,1) .* randn(1,length(Iinj));
            end

            Iinj(i,:)                  = Itmp;
            Iinj(i,(on./dt):(off./dt)) = Itmp(i,(on./dt):(off./dt)) + amplitude(i,1);

          else

            if test(1,2)==2
                delta = amplitude(1,2);
            end      

            Itmp = amplitude(1,1) .* randn(1,length(Iinj));
            Iinj(i,:)                  = Itmp;        
            Iinj(i,(on./dt):(off./dt)) = Itmp(1,(on./dt):(off./dt)) + amplitude(1,1) + (delta*(i-1));

          end


        end

      case 'sinusoid'
        t = (on./dt):(off./dt);
        tau = (30./dt); % ms
        Itmp = sin(2.*pi.*(t./tau));
        for i=1:numsweeps
            Iinj(i,(on./dt):(off./dt)) = Itmp + amplitude + (delta*(i-1));
        end

      case 'complexsin'
        t = (on./dt):(off./dt);
        components  = 2;
        frange      = 50;
        conversion  = 1000./dt;% From hertz to frequncy in dt
        for i=1:round(frange/components)
            amplitudes  = randn(1,1);
            offset      = randn(1,1);
            Itmp(i,:) = amplitudes .* (sin(2.*pi.*( ( [t./conversion] .* [(components-offset)*i] ) + (randn(1,1)*2*pi) ) ));
        end
        Iinj(1,(on./dt):(off./dt)) = amplitude*[sum(Itmp,1) ./ max(abs(sum(Itmp,1)))] + amplitude;
        if numsweeps>1
            for i=2:numsweeps
                Iinj(i,(on./dt):(off./dt)) = Iinj(1,(on./dt):(off./dt)) + (delta*(i-1));
            end
        end

      case 'broadband'

        [Itmp,t] = CreateBroadbandSignal(100,duration,dt);

        if test(1,1) > 1

            for i=1:numsweeps
                if test(1,2)==2
                    Iinj(i,(on./dt):(off./dt)) = (amplitude(i,2).*Itmp) + amplitude(i,1);
                else
                    Iinj(i,(on./dt):(off./dt)) = (amplitude(i,1).*Itmp) + amplitude(i,1);
                end
            end

        else

            for i=1:numsweeps
                if test(1,2)==2
                    Iinj(i,(on./dt):(off./dt)) = (amplitude(i,2).*Itmp) + amplitude(i,1) + (amplitude(i,2)*(i-1));
                else
                    Iinj(i,(on./dt):(off./dt)) = (amplitude(i,1).*Itmp) + amplitude(i,1);
                end
            end

        end

      otherwise
        error('Unrecognized stimulus type');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% NAME SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numsweeps>1
        basename = input('Please name this simulation:','s');
    else
        basename = 'single';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% EXECUTE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(Iinj(:,1));
%     [response,currents] = NeuronEngine(Vss,Iinj,dt,stochasticVec,conductanceVec,hdf,basename);
    [response,currents] = NeuronEngineFast(Vss,Iinj,dt,stochasticVec,conductanceVec,hdf,basename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SAVE GENERAL OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Save the call used for later reference.
    thecall = sprintf('Vss=%g,numsweeps=%g,duration=%g,on=%g,off=%g,IinjType=%s',Vss,numsweeps,duration,on,off,IinjType);
    thecurrents = cell(1,2);
    thecurrents{1,1}(1,:) = conductanceVec;
    thecurrents{1,1}(2,:) = stochasticVec;
    thecurrents{1,2}(:,:) = amplitude;

    save StoredDetails thecall thecurrents