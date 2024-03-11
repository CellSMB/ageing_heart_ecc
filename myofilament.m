function [VOI, STATES, ALGEBRAIC, CONSTANTS] = myofilament(tvec,ult_ca,maxt,xMHCa)
%   global model_t;global model_ca;
    global var_MHCa
    var_MHCa=xMHCa;
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel(tvec,ult_ca,maxt);
end

function [algebraicVariableCount] = getAlgebraicVariableCount()
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".
    algebraicVariableCount =35;
end
% There are a total of 4 entries in each of the rate and state variable arrays.
% There are a total of 16 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel(tvec,ult_ca,maxt)
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over;
    tspan = [0, maxt];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.01);
    
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS,tvec,ult_ca), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS,tvec,ult_ca);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI,tvec,ult_ca);

%     % Plot state variables against variable of integration
%     [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
%     figure();
%     plot(VOI, STATES);
%     xlabel(LEGEND_VOI);
%     l = legend(LEGEND_STATES);
%     set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('t in component environment (s)');
    
    %MYOFILAMENT FORCE MODULE
    LEGEND_CONSTANTS(:,1) = strpad('K_HTRPN_+ in component force (mMperms)');
    LEGEND_CONSTANTS(:,2) = strpad('K_HTRPN_- in component force (perms)');
    LEGEND_CONSTANTS(:,3) = strpad('K_LTRPN_+ in component force (mMperms)');
    LEGEND_CONSTANTS(:,4) = strpad('K_LTRPN_- in component force (perms)');
    
    LEGEND_CONSTANTS(:,5) = strpad('HTRPN_tot in component force (mM)');
    LEGEND_CONSTANTS(:,6) = strpad('LTRPN_tot in component force (mM)');
    
    LEGEND_CONSTANTS(:,7) = strpad('kpn_TRPN in component force (perms)');
    LEGEND_CONSTANTS(:,8) = strpad('SL_cort in component force (uM)');
    
    LEGEND_CONSTANTS(:,9) = strpad('fXB in component force (uM)');
    LEGEND_CONSTANTS(:,10) = strpad('gXB_min in component force (uM)');
    
    LEGEND_CONSTANTS(:,11) = strpad('zeta in component force (uM)');
    LEGEND_CONSTANTS(:,12) = strpad('gXB_min in component force (uM)');
    
    LEGEND_CONSTANTS(:,13) = strpad('VAM_max in component force (uM)');
    LEGEND_CONSTANTS(:,14) = strpad('KMAM_ATP in component force (uM)');
    LEGEND_CONSTANTS(:,15) = strpad('KiAM in component force (uM)');
    
    LEGEND_CONSTANTS(:,16) = strpad('Kpn_TRPN');
    
    LEGEND_ALGEBRAIC(:,16) = strpad('phi_myo');
    LEGEND_ALGEBRAIC(:,17) = strpad('f01');
    LEGEND_ALGEBRAIC(:,18) = strpad('f12');
    LEGEND_ALGEBRAIC(:,19) = strpad('f23');
    LEGEND_ALGEBRAIC(:,20) = strpad('g01');
    LEGEND_ALGEBRAIC(:,21) = strpad('g12');
    LEGEND_ALGEBRAIC(:,22) = strpad('g23');
    LEGEND_ALGEBRAIC(:,23) = strpad('g01_SL');
    LEGEND_ALGEBRAIC(:,24) = strpad('g12_SL');
    LEGEND_ALGEBRAIC(:,25) = strpad('g23_SL');
    
    LEGEND_ALGEBRAIC(:,26) = strpad('KCa_TRPN');
    LEGEND_ALGEBRAIC(:,27) = strpad('NTRPN');
    LEGEND_ALGEBRAIC(:,28) = strpad('Khalf_TRPN');
    LEGEND_ALGEBRAIC(:,29) = strpad('Knp_TRPN');
    
    LEGEND_ALGEBRAIC(:,30) = strpad('path_sum');
    LEGEND_ALGEBRAIC(:,31) = strpad('P1_max');
    LEGEND_ALGEBRAIC(:,32) = strpad('P2_max');
    LEGEND_ALGEBRAIC(:,33) = strpad('P3_max');
    
    LEGEND_STATES(:,1) = strpad('P0');
    LEGEND_STATES(:,2) = strpad('P1');
    LEGEND_STATES(:,3) = strpad('P2');
    LEGEND_STATES(:,4) = strpad('P3');
    LEGEND_STATES(:,5) = strpad('N0');
    LEGEND_STATES(:,6) = strpad('N1');
    LEGEND_STATES(:,7) = strpad('HTRPNCa');
    LEGEND_STATES(:,8) = strpad('LTRPNCa');
    
    LEGEND_STATES(:,9) = strpad('myo_Pim');
    LEGEND_ALGEBRAIC(:,14) = strpad('myo_cai');
    LEGEND_STATES(:,11) = strpad('myo_ATPi');
    LEGEND_STATES(:,12) = strpad('myo_ADP');
    
    LEGEND_ALGEBRAIC(:,34) = strpad('Norm Force');
    LEGEND_ALGEBRAIC(:,35) = strpad('Force');
    
    LEGEND_RATES(:,1) = strpad('dP0/dt');
    LEGEND_RATES(:,2) = strpad('dP1/dt');
    LEGEND_RATES(:,3) = strpad('dP2/dt');
    LEGEND_RATES(:,4) = strpad('dP3/dt');
    LEGEND_RATES(:,5) = strpad('dN0/dt');
    LEGEND_RATES(:,6) = strpad('dN1/dt');
    
    LEGEND_RATES(:,7) = strpad('dHTRPN/dt');
    LEGEND_RATES(:,8) = strpad('dLTRPN/dt');
    
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end
    
function [STATES, CONSTANTS] = initConsts()
    global var_MHCa; 

    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    
    CONSTANTS(:,1) = 100*(1e3);
    CONSTANTS(:,2) = (3.3e-4)*(1e3);
    CONSTANTS(:,3) = 100*(1e3);
    CONSTANTS(:,4) = (4e-2)*(1e3);
    CONSTANTS(:,5) = 0.14;
    CONSTANTS(:,6) = 0.7;
    
    CONSTANTS(:,7) = 0.04*(1e3);
    CONSTANTS(:,8) = 2.15;
    
    CONSTANTS(:,9) = 0.05*(1e3);
    CONSTANTS(:,10) = 0.1*(1e3);
    
    %percentage_alpha=(0.1*var_MHCa);
    %percentage_beta=0.9;
    %fa=0.2597;
    %ga=0.4639;
    %fb=0.0267;
    %gb=0.0596;
    
    %CONSTANTS(:,9) = (1e3)*(percentage_alpha*fa + percentage_beta*fb);
    %CONSTANTS(:,10) = (1e3)*(percentage_alpha*ga + percentage_beta*gb);
    
    CONSTANTS(:,11) = 0.1;
    CONSTANTS(:,12) = 0.1;
    
    CONSTANTS(:,13) = (7.2e-3)*(1e3);
    CONSTANTS(:,14) = 0.03;
    CONSTANTS(:,15) = 0.26;
    CONSTANTS(:,16) = 0.04*(1e3);
    
    STATES(:,1) = 1.6672012132596455e-3;
    STATES(:,2) = 1.4416741873793939e-3;
    STATES(:,3) = 2.6920336520753884e-3;
    STATES(:,4) = 2.3449372176401837e-3;
    STATES(:,5) = 9.9041768766118021e-1;
    STATES(:,6) = 1.4356068652011372e-3;
    
    STATES(:,7) = 1.3055735570840798e-1;
    STATES(:,8) = 1.8066410206828074e-2;
    
    %STATES(:,9) = 2.0;
    %STATES(:,10) = 0;
    %STATES(:,11) = 0;
    %STATES(:,12) = 0;
    
%     STATES=[0.705267474192312,0.346755423421977,0.00219128028633088,0.976768261761414];
    if (isempty(STATES)), warning('Initial values for states not set');, end
end    

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS,tvec,ult_ca)
    global algebraicVariableCount;%global model_ca; global model_t;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    
%     if abs(round(VOI) - VOI)
%         ALGEBRAIC(:,17) = (1e3)*ult_ca(1);
%     else
%         ALGEBRAIC(:,17) = (1e3)*interp1q(tvec,ult_ca,(rem(VOI,1000)*1e3+2e3));
%     end
    
    if abs(round(VOI) - VOI) == 0
        ALGEBRAIC(:,15) = ult_ca(1);
    else
        ALGEBRAIC(:,15) = interp1q(tvec,ult_ca,(rem(VOI*(1e3),1000)+2e3));
    end
    
    ALGEBRAIC(:,14) = ALGEBRAIC(:,15);
    
    %COMPUTE ALGEBRAICS AND RATES HERE
    
    ALGEBRAIC(:,16) = 1.000 + ((2.3000 - CONSTANTS(:,8))./power((2.3000 - 1.7000), 1.6000));
    
    ALGEBRAIC(:,17) = 3.*CONSTANTS(:,9);
    ALGEBRAIC(:,18) = 10.*CONSTANTS(:,9);
    ALGEBRAIC(:,19) = 7.*CONSTANTS(:,9);
    
    ALGEBRAIC(:,20) = CONSTANTS(:,10);
    ALGEBRAIC(:,21) = 2.*CONSTANTS(:,10);
    ALGEBRAIC(:,22) = 3.*CONSTANTS(:,10);
    
    ALGEBRAIC(:,23) = ALGEBRAIC(:,16).*CONSTANTS(:,10);
    ALGEBRAIC(:,24) = 2.*ALGEBRAIC(:,16).*CONSTANTS(:,10);
    ALGEBRAIC(:,25) = 3.*ALGEBRAIC(:,16).*CONSTANTS(:,10);
    
    ALGEBRAIC(:,26) = CONSTANTS(:,4)./CONSTANTS(:,3);
    ALGEBRAIC(:,27) = 3.5.*CONSTANTS(:,8) - 2.0000;
    ALGEBRAIC(:,28) = (1.000)./(1.000 + (ALGEBRAIC(:,26)./(1.4000e-3 - (0.8e-3).*((CONSTANTS(:,8) - 1.7000)./0.6000))));
    ALGEBRAIC(:,29) = CONSTANTS(:,16).*power(((STATES(:,8))./(ALGEBRAIC(:,28).*CONSTANTS(:,6))), ALGEBRAIC(:,27));
    
    ALGEBRAIC(:,30) = ALGEBRAIC(:,20).*ALGEBRAIC(:,21).*ALGEBRAIC(:,22) + ALGEBRAIC(:,17).*ALGEBRAIC(:,21).*ALGEBRAIC(:,22) + ALGEBRAIC(:,17).*ALGEBRAIC(:,18).*ALGEBRAIC(:,22) + ALGEBRAIC(:,17).*ALGEBRAIC(:,18).*ALGEBRAIC(:,19);
    ALGEBRAIC(:,31) = (ALGEBRAIC(:,17).*ALGEBRAIC(:,21).*ALGEBRAIC(:,22))./ALGEBRAIC(:,30);
    ALGEBRAIC(:,32) = (ALGEBRAIC(:,17).*ALGEBRAIC(:,18).*ALGEBRAIC(:,22))./ALGEBRAIC(:,30);
    ALGEBRAIC(:,33) = (ALGEBRAIC(:,17).*ALGEBRAIC(:,18).*ALGEBRAIC(:,19))./ALGEBRAIC(:,30);
    
    ALGEBRAIC(:,34) = (STATES(:,2) + STATES(:,6) + 2.000.*STATES(:,3) + 3.*STATES(:,4))./(ALGEBRAIC(:,31) + 2.*ALGEBRAIC(:,32) + 3.*ALGEBRAIC(:,33));
    ALGEBRAIC(:,35) = CONSTANTS(:,11).*ALGEBRAIC(:,34);
    
    %RATES
    
    RATES(:,1) = -(CONSTANTS(:,16) + ALGEBRAIC(:,17)).*STATES(:,1) + ALGEBRAIC(:,29).*STATES(:,5) + ALGEBRAIC(:,23).*STATES(:,2);
    
    RATES(:,2) = -(CONSTANTS(:,16) + ALGEBRAIC(:,18) + ALGEBRAIC(:,23)).*STATES(:,2) + ALGEBRAIC(:,29).*STATES(:,6) + ALGEBRAIC(:,17).*STATES(:,1) + ALGEBRAIC(:,24).*STATES(:,3);
    
    RATES(:,3) = -(ALGEBRAIC(:,19) + ALGEBRAIC(:,24)).*STATES(:,3) + ALGEBRAIC(:,18).*STATES(:,2) + ALGEBRAIC(:,22).*STATES(:,4);
    
    RATES(:,4) = -ALGEBRAIC(:,25).*STATES(:,4) + ALGEBRAIC(:,19).*STATES(:,3);
    
    RATES(:,6) = CONSTANTS(:,16).*STATES(:,2) - (ALGEBRAIC(:,29) + ALGEBRAIC(:,23)).*STATES(:,6);
    
    RATES(:,5) = -RATES(:,1) - RATES(:,2) - RATES(:,3) - RATES(:,4) - RATES(:,6);
    
    RATES(:,7) = ((CONSTANTS(:,1).*ALGEBRAIC(:,14).*(CONSTANTS(:,5) - STATES(:,7))) - (CONSTANTS(:,2).*STATES(:,7)));
    
    %RATES(:,7) = 1000*RATES(:,7);
    
    RATES(:,8) = (CONSTANTS(:,3).*ALGEBRAIC(:,14).*(CONSTANTS(:,6) - STATES(:,8))) - (CONSTANTS(:,4).*(1.0000 - ((2/3)*ALGEBRAIC(:,34))).*STATES(:,8));
    
    %RATES(:,8) = 1000*RATES(:,8);
    
    RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI,tvec,ult_ca)
% global model_ca;global model_t;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    
%     if abs(round(VOI) - VOI)
%         ALGEBRAIC(:,17) = (1e3)*ult_ca(1);
%     else
%         ALGEBRAIC(:,17) = (1e3)*interp1q(tvec,ult_ca,(rem(VOI,1000)*1e3+2e3));
%     end
    
    if abs(round(VOI) - VOI) == 0
        ALGEBRAIC(:,15) = ult_ca(1);
    else
        ALGEBRAIC(:,15) = interp1q(tvec,ult_ca,(rem(VOI*(1e3),1000)+2e3));
    end
    
    ALGEBRAIC(:,14) = ALGEBRAIC(:,15);
    
    %COMPUTE ALGEBRAICS HERE
    
    ALGEBRAIC(:,16) = 1.000 + ((2.3000 - CONSTANTS(:,8))./power((2.3000 - 1.7000), 1.6000));
    
    ALGEBRAIC(:,17) = 3.*CONSTANTS(:,9);
    ALGEBRAIC(:,18) = 10.*CONSTANTS(:,9);
    ALGEBRAIC(:,19) = 7.*CONSTANTS(:,9);
    
    ALGEBRAIC(:,20) = CONSTANTS(:,10);
    ALGEBRAIC(:,21) = 2.*CONSTANTS(:,10);
    ALGEBRAIC(:,22) = 3.*CONSTANTS(:,10);
    
    ALGEBRAIC(:,23) = ALGEBRAIC(:,16).*CONSTANTS(:,10);
    ALGEBRAIC(:,24) = 2.*ALGEBRAIC(:,16).*CONSTANTS(:,10);
    ALGEBRAIC(:,25) = 3.*ALGEBRAIC(:,16).*CONSTANTS(:,10);
    
    ALGEBRAIC(:,26) = CONSTANTS(:,4)./CONSTANTS(:,3);
    ALGEBRAIC(:,27) = 3.5.*CONSTANTS(:,8) - 2.0000;
    ALGEBRAIC(:,28) = (1.000)./(1.000 + (ALGEBRAIC(:,26)./(1.4000e-3 - 0.8e-3.*((CONSTANTS(:,8) - 1.7000)./0.6000))));
    ALGEBRAIC(:,29) = CONSTANTS(:,16).*power(((STATES(:,8))./(ALGEBRAIC(:,28).*CONSTANTS(:,6))), ALGEBRAIC(:,27));
    
    ALGEBRAIC(:,30) = ALGEBRAIC(:,20).*ALGEBRAIC(:,21).*ALGEBRAIC(:,22) + ALGEBRAIC(:,17).*ALGEBRAIC(:,21).*ALGEBRAIC(:,22) + ALGEBRAIC(:,17).*ALGEBRAIC(:,18).*ALGEBRAIC(:,22) + ALGEBRAIC(:,17).*ALGEBRAIC(:,18).*ALGEBRAIC(:,19);
    ALGEBRAIC(:,31) = (ALGEBRAIC(:,17).*ALGEBRAIC(:,21).*ALGEBRAIC(:,22))./ALGEBRAIC(:,30);
    ALGEBRAIC(:,32) = (ALGEBRAIC(:,17).*ALGEBRAIC(:,18).*ALGEBRAIC(:,22))./ALGEBRAIC(:,30);
    ALGEBRAIC(:,33) = (ALGEBRAIC(:,17).*ALGEBRAIC(:,18).*ALGEBRAIC(:,19))./ALGEBRAIC(:,30);
    
    ALGEBRAIC(:,34) = (STATES(:,2) + STATES(:,6) + 2.000.*STATES(:,3) + 3.*STATES(:,4))./(ALGEBRAIC(:,31) + 2.*ALGEBRAIC(:,32) + 3.*ALGEBRAIC(:,33));
    ALGEBRAIC(:,35) = CONSTANTS(:,11).*ALGEBRAIC(:,34);
    
    
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end

    
    
    
    
    
    
    