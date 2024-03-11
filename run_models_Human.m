function [out, tvec, ult_ca, tvec_F, ult_F]=run_models_Human(params)

%params = kp;
[VOI, STATES, ALGEBRAIC, CONSTANTS] = IMW_Human_ECC_wIP3(params,[]);
    
period=CONSTANTS(:,12); %perdio in ms; default is 1000
% stability
numOsc=round(VOI(end)/period);
penUlt=interp1(VOI-(numOsc-2)*period,STATES(:,3)*1e3,0:(period-1));
ult=interp1(VOI-(numOsc-1)*period,STATES(:,3)*1e3,0:(period-1));
oscDiff=sum((ult-penUlt).^2/range(ult).^2);
redoNum=0;
maxRedo=50;
while oscDiff>1e-3 && redoNum<maxRedo
        
    [VOI, STATES, ALGEBRAIC, CONSTANTS] = IMW_Human_ECC_wIP3(params,STATES(end,:));
        
     redoNum=redoNum+1;
     penUlt=interp1(VOI-(numOsc-2)*period,STATES(:,3)*1e3,0:(period-1));
     ult=interp1(VOI-(numOsc-1)*period,STATES(:,3)*1e3,0:(period-1));
     oscDiff=sum((ult-penUlt).^2/range(ult).^2);
end
if redoNum>=maxRedo
    warning(strcat('Not great clarity'));
end

%end stability
lastTransient=find(VOI>=(VOI(end)-period));
out = getQualsHuman(VOI, STATES, ALGEBRAIC, CONSTANTS,lastTransient);

ultimate_transient = VOI>(VOI(end)-CONSTANTS(:,12));
penultimate_transient = VOI>(VOI(end)-2*CONSTANTS(:,12))&VOI<=(VOI(end)-CONSTANTS(:,12));
ultimate_transient_cai = STATES(ultimate_transient,3);
penultimate_transient_cai = STATES(penultimate_transient,3);
transient_time=VOI(ultimate_transient);

%assign outputs
tvec=transient_time;
ult_ca=ultimate_transient_cai;

%FORCE MODULE
maxt = 100;
xMHCa=params(10);
[VOI, STATES, ALGEBRAIC, CONSTANTS] = myofilament(tvec,ult_ca,maxt,xMHCa);
%out3=[STATES(end,:)];
%alg_myo = ALGEBRAIC;
%states = STATES;

period=1;
lastTwitch=find(VOI>=(VOI(end)-period)); %get the last twitch 
out3=getForceQuals(VOI, STATES, ALGEBRAIC, CONSTANTS,params,lastTwitch); %find the twitch qualities

ultimate_transient = VOI>(VOI(end)-period);
penultimate_transient_F = VOI>(VOI(end)-2*period)&VOI<=(VOI(end)-period);
ultimate_transient_F = ALGEBRAIC(ultimate_transient,35);
penultimate_transient_F = ALGEBRAIC(penultimate_transient_F,35);
transient_time_F=VOI(ultimate_transient);

tvec_F=transient_time_F;
ult_F=ultimate_transient_F;

out=[out out3];

% maxt = 6e3;
% [~, STATES, ~, ~] = NFAT_Cooling2009(tvec,ult_ca,maxt);
% out2=[STATES(end,:)];
% out=[out out2];

%maxt = 3.6e3; %3.6;
%num_reps=floor(maxt*1000/period)+1; %floor(maxt*1000/period)+1;