%% Generate caTs for range of parameters
clear all; close all; clc;

N_SERCA=4;
N_MHC=1;% total no. of simulated model output to quantify sensitivity
kp=[0 0 1 1 1 1 1 1 1 1]; %0 0 1 1 0.6680 1 1 1 1]; %; 0.1 0.1 1 1 1 1 1 1 1; 0 0 1 1 0.8 1 1 1 1; 0.1 0.1 1 1 0.8 1 1 1 1];
n_p=size(kp,2);
%PS=[];
SERCA_vars=linspace(1,0.668,N_SERCA); %0.668
MHC_vars=linspace(1,0.7,N_MHC);
PS=[];
output_quals=[];
tvec_all={};
ult_ca_all={};
tvec_F_all={};
ult_F_all={};

%for i = 1:n_p
    %PS=[PS ((1e-5 + 2.*rand(N,1)))];
    %PS=[PS ((1 + (1e-5).*rand(N,1)))];
%end

%testing parameters set
%PS=[0 0 1 1 1 1 1 1 1];

for i=1:N_SERCA
    i
    %reset variables
    kp=[0 0 1 1 1 1 1 1 1 1]
    
    kp(5)=SERCA_vars(i)
    
    for j=1:N_MHC
    
        kp(10)=MHC_vars(j)
        
        [out, tvec, ult_ca, tvec_F, ult_F]=run_models_Human(kp);
        PS=[PS;kp];
        output_quals=[output_quals; out];
        tvec_all{end+1}=tvec;
        ult_ca_all{end+1}=ult_ca;
        tvec_F_all{end+1}=tvec_F;
        ult_F_all{end+1}=ult_F;
    end
end