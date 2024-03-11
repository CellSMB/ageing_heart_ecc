function results=getForceQuals(VOI, STATES, ALGEBRAIC, CONSTANTS,params,lastTransient)
   period = 1;
   
   %dF/dt
   dFdt_max = max(diff(ALGEBRAIC(lastTransient,35))./diff(VOI(lastTransient)));
   dFdt_min = min(diff(ALGEBRAIC(lastTransient,35))./diff(VOI(lastTransient)));
    
   %average calcium
   avCa = mean(ALGEBRAIC(lastTransient,35));
    
   % baseline
   baseline=min(ALGEBRAIC(lastTransient,35));
   %concbaselineF=min(STATES(lastTransient,1));
   %transientFluo4 = STATES(lastTransient,1)/concbaselineF;
   
   % amplitude & ttpeak
   [amplitude,peakPos]=max(ALGEBRAIC(lastTransient,35));
   ttpeak = VOI(lastTransient(1)+peakPos-1)-VOI(lastTransient(1));
   %[amplitudeF,peakPosF]=max(transientFluo4);
   %ttpeakF = VOI(lastTransient(1)+peakPosF-1)-VOI(lastTransient(1));
   
   % Duration 50%
   threshold = (amplitude+baseline)/2;
   %thresholdF = (amplitudeF+1)/2;
   dur50pos = find(ALGEBRAIC(lastTransient,35)>=threshold);
   dur50 = VOI(dur50pos(end)+lastTransient(1)-1)-VOI(dur50pos(1)+lastTransient(1)-1);
   %dur50Fpos = find(transientFluo4>=thresholdF);
   %dur50F = VOI(dur50Fpos(end)+lastTransient(1)-1)-VOI(dur50Fpos(1)+lastTransient(1)-1);
  
   % Duration 90%
   threshold = 0.1*amplitude+0.9*baseline;
   %thresholdF = amplitudeF*0.1+0.9;
   dur90pos = find(ALGEBRAIC(lastTransient,35)>=threshold);
   dur90 = VOI(dur90pos(end)+lastTransient(1)-1)-VOI(dur90pos(1)+lastTransient(1)-1);
   %dur90Fpos = find(transientFluo4>=thresholdF);
   %dur90F = VOI(dur90Fpos(end)+lastTransient(1)-1)-VOI(dur90Fpos(1)+lastTransient(1)-1);
   
   % Duty cycle
   dutyCycle = trapz(VOI(lastTransient),ALGEBRAIC(lastTransient,35)-baseline)/((amplitude-baseline)*period);
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
   results = [avCa,amplitude,1000*dur50,1000*dur90,1000*ttpeak,baseline,dutyCycle,dFdt_max,dFdt_min];
end