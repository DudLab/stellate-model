function [Vss] = LookupSSvoltage(conductanceVec,Iinj);
% Get the ss voltage from a given injected current for a model with
% conductances set by the conductance Vec

% Create the net current table
[currents,voltage] = GetSSIV(-100.01342567846234212890,-40.01342567846234212890,0.025,conductanceVec);

% Find the net current that matches Iinj
if (Iinj==0)
    voltageindex = find(currents(9,:)>=0);
else
    voltageindex = find(currents(9,:)>=Iinj);
end

if length(voltageindex)==0
    Vss = -52.1
    disp('You have injected too much current to estimate steady-state, using default of -52.1');
else
    Vss = voltage(1,voltageindex(1,1))
end

% Retrieve the corresponding voltage
