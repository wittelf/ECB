%% Empirical charge buildup (ECB)model
% Falk Wittel and Eashan Saikia, ETH Zurich May 2020
% Evaluates the charge buildup in the flytrap for single & multiple stimuli
% ---------------------------------------------------------------------------
% Inputs     : Mechanical Stimulus in terms of Angular Deflection and Angular Velocity
% Outputs    : Number of Action Potential (AP) and predicts trap closure.
% ----------------------------------------------------------------------
% Model Parameters:-
%   ----------------------
%   Name              Units        Description
dt      = 0.00005;	  %   s          : Fixed time increment (delta t)
tau_RP	= 0.640460;   %   s          : Decay constant for receptor potential (RP)
tau_STM	= 39.713477;  %   s          : Decay constant for Short term memory (STM)
k       = 40.203466;  %   uC/rad     : Sensitivity parameter or mechanotransduction
a=1.994280;b=3.477958;%    -          : Desensitivity parameters for (k) after an AP
Qth_RP  = 1;          %   uC         : Threshold charge for RP to elicit an AP
Qth_STM = 10.0;       %   uC         : Threshold charge for STM to initiate trap closure
QAP     = 6;          %	  uC         : Fixed charge added into STM after an AP
maxflux = 0.004651;	  %	  mC/sec     : Limiting charge flux
t_RP    = 1.029194;   %	  s          : Refractory period during which no AP is possible
RF      = 1;          %	-          : Factor for considering desensitization after first firing

%   Simulation parameter    
poke    = 3;        %	-          : No. of stimuli (deflections)
t_gap   = 1;        %	s          : Wait time between 2 consecutive pokes
w       = 1;        %   rad/s      : Set angular velocity here example 1
phi     = 0.2;      %   rad        : Set angular deflection here example 0.2

%   Variables:
T       = 2*phi/w;	%	s           : Time period of one stimulus (advance+return)
triggered = false;  %   -           : false if trap remains open, true if closed
w_vec   = w*ones(1,round(T/dt));        %:  angular velocity vector
w_vec   = [w_vec, zeros(1,round(t_RP/dt))]; %Add refractory period after it
for i = 2:poke; w_vec = [w_vec, zeros(1,round(t_gap/dt)), w*ones(1,round(T/dt))]; end; 
w_vec = [w_vec, zeros(1,round(5*t_RP/dt))];%5 times waiting afterwards
Q_RP     = zeros(size(w_vec)); %uC	: Charge in the RP
Q_STM    = zeros(size(w_vec)); %uC	: Charge in the STM
nap      = 0;       %   -           : Number of triggered APs
timenew  = t_RP;	%   s           : time since last firing 

%   Functions:
mf      = @(x) 1./(1 + exp(-5.*(x)));                        % function handle for sigmoid function to control maximal flux
mrm     = @(C,a,b) 1- sigmf(C,[a b]);% function handle for desensitization after each stimulus (sigmoidal membership function)
%  Main loop
for m = 1:length(w_vec) - 1     %loop over all the time steps 
    timenew = timenew + dt;            %increment time since last firing
    added_charge = k * dt * w_vec(m) * mrm(Q_RP(m),a,b);                       %calculate the incremental charge for the RP
    added_charge = min([mf(added_charge/maxflux)*maxflux, added_charge]);   %restriction for max.flux
    Q_RP(m+1) = Q_RP(m) * exp(-dt/tau_RP) + added_charge;   %decay and increase of RP 
    Q_STM(m+1) = Q_STM(m) * exp(-dt/tau_STM);               %decay of STM 
    
    %check if AP firing threshold is reached
    if (Q_RP(m+1) >= Qth_RP) && (timenew > t_RP)         
        nap = nap + 1;                                % increment AP counter
        timenew = 0;                                  % reset time since last firing
        myRF = 1;                                     %Proviede correct 
        if (nap >= 1);    myRF = RF;  end       %After 1st AP. Could be in the 1st or at a later poke
        Q_STM(m+1) = Q_STM(m+1)  + RF * QAP;    %Add reduced AP charge to the STM                  
    end
    %check for trap closure  
    if (Q_STM(m+1) >= Qth_STM);   triggered = 1; end        
end

display(['Trap closed =',num2str(triggered),', ',num2str(nap),' APs triggered, Final RP charge: ',num2str(Q_RP(end)),...
    ' Final STM charge: ',num2str(Q_STM(end)),' Max. STM charge: ',num2str(max(Q_STM))])