##########################################################################################################
#
#   <This program implements the Dynamic Computational Model of Motivation (DCMM)
#   which is a conectionist model of motivation based on SDT and HMIEM. The framework
#   consists on interconected Continuos Attractor Neural Netoworks (CANN)>
#
#   Copyright (C) <2018>  <Hendry Ferreira Chame>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, version 3.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#    
#   Programmed by: Hendry Ferreira Chame
#   Institution: Federal University of Rio Grande (FURG), Rio Grande, BRAZIL
#   E-mail: hendrychame@furg.br, hendryf@gmail.com
#     
#   Paper:  A dynamic computational model of motivation based on self-determination theory and CANN
#   Authors: Hendry Ferreira Chame, Fernanda Pinto Mota, Silvia Silva da Costa Botelho
#   Information Sciences, Elsevier, 2018. https://doi.org/10.1016/j.ins.2018.09.055
#
#   Description:  this program implenets Case Study 2 (CS2). 
#                 The plots generated are in the folder 'plots' 
#
##########################################################################################################


clc;
close all;
clear all;

display("#################################################################################################");
display("")
display("  <GNU Octave Implementation of DCMM>  Copyright (C) <2018>  <Hendry Ferreira Chame>          ");
display("  This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome    "); 
display("  to redistribute it under certain conditions; for details consult the file 'LICENSE.txt'.      ");
display("")
display("  For more details about the framework, please consult the following reference:                 "); 
display("")
display(" \"A dynamic computational model of motivation based on self-determination theory and CANN\"    ");  
display("  Hendry Ferreira Chame, Fernanda Pinto Mota, Silvia Silva da Costa Botelho                     ");
display("  Information Sciences, Elsevier, 2018                                                          ");
display("")
display("#################################################################################################");
display("")


display ("DCMM program stared!")
display ("Begin of Case study 3 (CS3) - The field study")
fflush(stdout); 


%% Shared simulation  parameters

motNeurons = 6;    % Number motivation layer units
nContext = 4;      % Life context simulated (1: Education, 2: Interpersonal relations, 3: Leissure, 4: Null)
nMediators = 3;    % Number of mediator processes (1: Autonomy, 2: Competence, 3: Relatedness)
medNeurons = 5;    % Number of mediator layers units
nConsequences = 3; % Number of consequences outputs (1: Behavior, 2: Cognition, 3: Affection)
conNeurons = 10;   % Number of consequence layer units
  
sigmaCann = 0.3;   % Sigma CANN gaussian activation
sigmaMed = 0.45;   % Sigma Mediators gaussian activation
sigmaCon = 1.0;    % Sigma Consequences gaussian activation

%% DCMM parameter setting 

setWeightParameters(sigmaCann, sigmaMed, motNeurons, medNeurons);
 
load("parameters/WAuToM.m");
load("parameters/WCoToM.m");
load("parameters/WReToM.m");
load("parameters/WCANN.m");
load "measures/moGlobal.m";
load "measures/moContextual.m";
load "measures/moSituational.m";
load "measures/situationalNeeds.m";
load "measures/situationalNeedsImp.m";
load "measures/selfReport.m";

display("Computing CS3 ... " );
fflush(stdout); 

dt = 5.0; # in sec

% Parameters for situational level 
tS = dt;
eS = 0.1;
g1 = g2 = g3 = 1.0;
g4 = 0.5;

% Parameters for contextual level 
tC = tS;  
eC = eye(nContext).*0.1;
g5 = g6 = g7 = g8 = g9 = gAMS = zeros(nContext,1);
g9(1) = 0.5;
gAMS(1) = 2;

% Parameters for global level 
tG = tC;
eG = 0.1;
g10 = g11 = g12 = g13 = 0;
gGMS = 10;

%% Inserting information to DCMM from the scales
 
n = size(selfReport, 1); % number of subjects

nGMSstates = 5; % number of GMS motivation states
nAMSstates = 5; % number of AMS motivation states
nSIMSstates = 4; % number of SIMS motivation states
nMedFactors = nMediators*2 ; % satisfaction and frustration 

% database statistics index-> (1:Mean,2:Median,3:Mode,4:Max,5:Min,6:stDev)
index = 1; % mean value in scales subcategories is used
  
% getting data according to the statistics used 

i1 = (index-1)*nGMSstates + 1;
i2 = (index-1)*nAMSstates + 1;
i3 = (index-1)*nSIMSstates + 1;

scale_GMS = moGlobal(:,i1:i1+nGMSstates-1);
scale_AMS = moContextual(:,i1:i1+nAMSstates-1);
scale_SMS = moSituational(:,i2:i2+nSIMSstates-1);
scale_BPNSF = situationalNeeds(:,i3:i3+nMedFactors-1);
scale_BPNSF_i = situationalNeedsImp(:,i3:i3+nMedFactors-1);

  
% Encapsulating Motivation scales into motivation layers by linear projection

h_GMS = [scale_GMS(:,1:4) mean(scale_GMS(:,4:5),2) scale_GMS(:,5)]';
h_GMS = h_GMS ./ sum(h_GMS,1); %normalizing 
    
h_AMS = [scale_AMS(:,1:4) mean(scale_AMS(:,4:5),2) scale_AMS(:,5)]';
h_AMS = h_AMS ./ sum(h_AMS,1); %normalizing

h_SMS = [scale_SMS(:,1:2) mean(scale_SMS(:,2:3),2) scale_SMS(:,3) mean(scale_SMS(:,3:4),2) scale_SMS(:,4)]';    
h_SMS = h_SMS ./ sum(h_SMS,1); %normalizing

%% --- situational Mediators
      
% Loading BPNSF rates in mediators layers

% It is calculated as the mean between the rates obtained in the items measuring satisfaction,
% plus the inversed rates obtained in the items measuring frustration.

%% setting Mediator's activation
step =  (motNeurons-1)/(medNeurons-1);
preferedUnitValue = [1:step:motNeurons]';  % preferred value encoded by the units
  
actFunForNStepsA = mean([scale_BPNSF(:,1) (medNeurons+1)-scale_BPNSF(:,2)],2)';       
actFunForNStepsC = mean([scale_BPNSF(:,3) (medNeurons+1)-scale_BPNSF(:,4)],2)';
actFunForNStepsR = mean([scale_BPNSF(:,5) (medNeurons+1)-scale_BPNSF(:,6)],2)';

f_As = generateGaussian(sigmaMed,actFunForNStepsA,preferedUnitValue);
f_Cs = generateGaussian(sigmaMed,actFunForNStepsC,preferedUnitValue);
f_Rs = generateGaussian(sigmaMed,actFunForNStepsR,preferedUnitValue);

%% --- Contextual Mediators

f_Ac_k = zeros(medNeurons,n,nContext);
f_Cc_k = zeros(medNeurons,n,nContext);
f_Rc_k = zeros(medNeurons,n,nContext); 

%% --- Global Mediators

f_Ag = f_Cg = f_Rg = zeros(medNeurons,n);

%% Setting initial states

% Activation Ms
h_Ms = ones(motNeurons,n).*(1/motNeurons); % uniform distribution
h_Ms_all = [];

% Activation Mc
h_Mc = ones(motNeurons,n).*(1/motNeurons); % uniform distribution
h_Mc_all = [];

% Activation Mc for the k contexts
h_Mc_K = ones(motNeurons,n,nContext).*(1/motNeurons); % uniform distribution
h_Mc_K_all =  [];
  
% Activation Mg  
h_Mg = ones(motNeurons,n).*(1/motNeurons); % uniform distribution
h_Mg_all = [];

% Applying the inhibition factors
W_Ms = WCANN .+ eS; 
W_Mg = WCANN .+ eG; 

% inter-contexts recurrent connexions 
W_Mc_K = zeros(motNeurons*nContext,motNeurons,nContext);
index = [0: nContext].*motNeurons;

for i = 2 : nContext+1
  W_Mc_K(index(i-1)+1:index(i),1:motNeurons,i-1) = WCANN;  
endfor 
  
% CAL connections and activation 
WCALtoMc_K = zeros(motNeurons,nContext,nContext);
% unbiased contribution from contexts
WCALtoMc_K(:,1,1) = ones(motNeurons,1);
WCALtoMc_K(:,2,2) = ones(motNeurons,1);
WCALtoMc_K(:,3,3) = ones(motNeurons,1);
WCALtoMc_K(:,4,4) = ones(motNeurons,1);

% Only education is relevant for the study
h_CAL = [1 0 0 0]';
h_CAL = h_CAL./sum(h_CAL);

% TNL connections and activation
WTNLtoAs = WTNLtoCs = WTNLtoRs = zeros(medNeurons,nMediators);
% unbiased contribution from mediators
WTNLtoAs(:,1) =  ones(medNeurons,1);
WTNLtoCs(:,2) =  ones(medNeurons,1);
WTNLtoRs(:,3) =  ones(medNeurons,1);

f_TNL = zeros(nMediators,n); 
f_TNL(1,:) = mean([scale_BPNSF_i(:,1) (medNeurons+1)-scale_BPNSF_i(:,2)],2);       
f_TNL(2,:) = mean([scale_BPNSF_i(:,3) (medNeurons+1)-scale_BPNSF_i(:,4)],2);       
f_TNL(3,:) = mean([scale_BPNSF_i(:,5) (medNeurons+1)-scale_BPNSF_i(:,6)],2);       
f_TNL = f_TNL ./ sum(f_TNL,1);

%%%%%%%%%

%% Prediction begin 

% Global level

f_Mg  = h_Mg;
   
h_Mg = (1-(dt/tG))*h_Mg + (dt/tG)*(W_Mg*f_Mg + g10*WAuToM*f_Ag + g11*WCoToM*f_Cg + g12*WReToM*f_Rg + g13*h_Mc + + gGMS*h_GMS);

h_Mg = h_Mg./sum(h_Mg,1);
h_Mg_all = [h_Mg_all h_Mg];
  
% DCMM computation

h_Mc = zeros(motNeurons,n); 

for k = 1 : nContext
        
  W_Mc_k = reshape (W_Mc_K(:,:,k), motNeurons*nContext,motNeurons);
  f_Mc_k = zeros(motNeurons,n);
  
  for y = 1 : nContext
  
    iL = (y-1)*motNeurons+1;
    rL = iL + motNeurons -1;
    W_Mc_y_k = W_Mc_k(iL:rL,:) .+ eC(k,y);
    f_Mc_y =  h_Mc_K(:,y);
    f_Mc_k  = f_Mc_k + W_Mc_y_k*f_Mc_y;
    
  endfor

  f_Ac = reshape (f_Ac_k(:,:,k), medNeurons,n);
  f_Cc = reshape (f_Cc_k(:,:,k), medNeurons,n);
  f_Rc = reshape (f_Rc_k(:,:,k), medNeurons,n);
  
  h_Ac  = (WAuToM*f_Ac);
  h_Cc  = (WCoToM*f_Cc);
  h_Rc  = (WReToM*f_Rc);
  
  h_Mc_k = (1-(dt/tC))*h_Mc_K(:,k) + (dt/tC)*(f_Mc_k + g5(k)*h_Ac + g6(k)*h_Cc + g7(k)*h_Rc + g8(k)*h_CAL(k)*h_Ms + g9(k)*h_CAL(k)*h_Mg + gAMS(k)*h_AMS);
  
  h_Mc_k = h_Mc_k./(sum(h_Mc_k,1)+1.0e-10);   
  h_Mc_K(:,:,k) = h_Mc_k;
  
  WCALtoMc_k = reshape (WCALtoMc_K(:,:,k), motNeurons,nContext);     
  h_Mc = h_Mc + WCALtoMc_k*h_CAL.*h_Mc_k;  
  
endfor

h_Mc_K_all = [h_Mc_K_all h_Mc_K];
h_Mc  = h_Mc./(sum(h_Mc));
h_Mc_all = [h_Mc_all h_Mc];
  
% Situational level

f_Ms  = h_Ms;

h_As  = WTNLtoAs*f_TNL.*f_As;
h_Cs  = WTNLtoCs*f_TNL.*f_Cs;
h_Rs  = WTNLtoRs*f_TNL.*f_Rs;
   
h_Ms = (1-(dt/tS))*h_Ms + (dt/tS)*(W_Ms*f_Ms + g1*WAuToM*h_As + g2*WCoToM*h_Cs + g3*WReToM*h_Rs + g4*h_Mc);

h_Ms = h_Ms./sum(h_Ms,1);
h_Ms_all = [h_Ms_all h_Ms];
   
%% Predicting consequences
   
step = (conNeurons-1)/(motNeurons-1);
preferedUnitValue = [1:step:conNeurons]'; % preferred value encoded by the units

%% Weight matrix for linear projection from motivation-like layer to consequences layer

WMtoCons = zeros(conNeurons,motNeurons);
consScale = [1:1:conNeurons];

for i=1:motNeurons
  WMtoCons(:,i) = generateGaussian(sigmaCon,preferedUnitValue(i),consScale);
end
WMtoCons = WMtoCons./ sum(WMtoCons,1);

%% Predictions

predDCMM = WMtoCons*h_Ms;
predSIMS = WMtoCons*h_SMS;
predAMS  = WMtoCons*h_AMS;
predGMS  = WMtoCons*h_GMS;

%% select he most activated unit by the winner-takes-all policy
[idx wtaDCMM] = max(predDCMM,[],1); % max
[idx wtaSIMS] = max(predSIMS,[],1); % max 
[idx wtaAMS]  = max(predAMS,[],1); % max 
[idx wtaGMS]  = max(predGMS,[],1); % max
  
   
%% calculating the Pearson product-moment correlation coefficient between prediction and subjects' self-report
pearson = [corr(wtaDCMM,selfReport) ; corr(wtaSIMS,selfReport) ; corr(wtaAMS,selfReport); corr(wtaGMS,selfReport)];

display("Pearson product-moment correlation coefficient:");
display("P/C \t[Cognition - Behavior - Affect]");
display("----------------------------------------------");
display(strcat("DCMM \t ", num2str(pearson(1,1),"%1.4f"), "\t", num2str(pearson(1,2),"%1.4f"), "\t", num2str(pearson(1,3),"%1.4f")));
display(strcat("SIMS \t ", num2str(pearson(2,1),"%1.4f"), "\t", num2str(pearson(2,2),"%1.4f"), "\t", num2str(pearson(2,3),"%1.4f")));
display(strcat("AMS \t ", num2str(pearson(3,1),"%1.4f"), "\t", num2str(pearson(3,2),"%1.4f"), "\t", num2str(pearson(3,3),"%1.4f")));
display(strcat("GMS \t ", num2str(pearson(4,1),"%1.4f"), "\t", num2str(pearson(4,2),"%1.4f"), "\t", num2str(pearson(4,3),"%1.4f")));
display("----------------------------------------------");

%% Plotting the obtained correlations 

figure(1);
s = size(pearson,2);
plot([1:s], pearson(1,:), 'b');
hold on;
plot([1:s], pearson(2,:), 'r');
plot([1:s], pearson(3,:), 'k');
plot([1:s], pearson(4,:), 'g');
%hx = xlabel("Output");
hy = ylabel("Correlation");
ht = title("Pearson product-moment correlation coefficient");  
legend("DCMM ","SIMS ","AMS ","GMS ", "location", "northwest");
countries = {'Cognition','Behavior','Affect'};
set(gca,'XTick',1:3,'XTickLabel',countries)
print -deps -color plots/CS2_pearson.eps

  
display ("End of Case study 3 (CS3) - The field study");
display ("DCMM program finished")

 
