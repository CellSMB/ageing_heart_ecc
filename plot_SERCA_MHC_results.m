%%
%Correlation coefficients
PS=[PS,log2(PS(:,5))];

ccs_all=corrcoef([PS(:,11),output_quals]);
ccs_interest=ccs_all(1,[9,15,16]);

%peak, baseline

%white background, black and red graphs 

%%
N_SERCA=4;
SERCA_range=[min(log2(SERCA_vars)),max(log2(SERCA_vars))];
N_MHC=1;
MHC_range=[min(log2(MHC_vars)),max(log2(MHC_vars))];

figure
hold on
for i=1:N_SERCA*N_MHC
    
    if 1
        
        %SERCA plot
        
        subplot(1,2,1)
        hold on
        plot(tvec_all{i}-2000,1000*ult_ca_all{i});
    
        
        if i==1
            xlabel('time(ms)')
            ylabel('Intracellular [Ca^{2+}]')
        end
        
        %Force plot
        
        subplot(1,2,2)
        hold on
        plot(1000*(tvec_F_all{i}-99),ult_F_all{i});
    
        
        if i==1
            xlabel('time(ms)')
            ylabel('Force (A.U.)')
        end
               
    end

end

