function results=getQualsHuman(VOI, STATES, ALGEBRAIC, CONSTANTS,lastTransient)
   period = 1;
    
   %average calcium
   avCa = mean(STATES(lastTransient,3));
    
   % baseline
   baseline=min(STATES(lastTransient,3));
   %concbaselineF=min(STATES(lastTransient,1));
   %transientFluo4 = STATES(lastTransient,1)/concbaselineF;
   
   % amplitude & ttpeak
   [amplitude,peakPos]=max(STATES(lastTransient,3));
   ttpeak = VOI(lastTransient(1)+peakPos-1)-VOI(lastTransient(1));
   %[amplitudeF,peakPosF]=max(transientFluo4);
   %ttpeakF = VOI(lastTransient(1)+peakPosF-1)-VOI(lastTransient(1));
   
   % Duration 50%
   threshold = (amplitude+baseline)/2;
   %thresholdF = (amplitudeF+1)/2;
   dur50pos = find(STATES(lastTransient,3)>=threshold);
   dur50 = VOI(dur50pos(end)+lastTransient(1)-1)-VOI(dur50pos(1)+lastTransient(1)-1);
   %dur50Fpos = find(transientFluo4>=thresholdF);
   %dur50F = VOI(dur50Fpos(end)+lastTransient(1)-1)-VOI(dur50Fpos(1)+lastTransient(1)-1);
  
   % Duration 90%
   threshold = 0.1*amplitude+0.9*baseline;
   %thresholdF = amplitudeF*0.1+0.9;
   dur90pos = find(STATES(lastTransient,3)>=threshold);
   dur90 = VOI(dur90pos(end)+lastTransient(1)-1)-VOI(dur90pos(1)+lastTransient(1)-1);
   %dur90Fpos = find(transientFluo4>=thresholdF);
   %dur90F = VOI(dur90Fpos(end)+lastTransient(1)-1)-VOI(dur90Fpos(1)+lastTransient(1)-1);
   
   % Duty cycle
   dutyCycle = trapz(VOI(lastTransient),STATES(lastTransient,3)-baseline)/((amplitude-baseline)*period);
   %dutyCycleF = trapz(VOI(lastTransient),transientFluo4-1)/((amplitudeF-1)*period);
  
   %clear x_sub x y
   %figure
   %plot(x,y_norm,'b')
   %hold on
   %expfcn = @(x) exp(x./-B(1));
   %plot(x,expfcn(x),'r-')
   
   % Init states
   inits=STATES(end,:);
   %Fluxes
   %I_IP3R
   %maxip3r=max(ALGEBRAIC(lastTransient,66));
   %maxryr=max(ALGEBRAIC(lastTransient,35));
   results = [avCa,amplitude,dur50,dur90,ttpeak,baseline,dutyCycle];
end