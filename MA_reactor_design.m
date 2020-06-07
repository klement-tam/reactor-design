%% README
%The way this code is written is that instead of plotting everything at once (i.e. Adiabatic, Isothermal, sensitivity analysis for
% b1, b2, n etc.), it uses dialogue prompts to allow the user to select first the type of system, and second which parameter is desired
% for sensitivity analysis, thereafter plotting the relevant graphs. This gives the code more flexibility and efficiency.

%If sensitivity analysis for n-butane:air ratio, Inlet Temperature or Inlet Pressure is selected, an additional dialogue box will 
% pop up asking the user if Further Optimisation is desired. These are inlet parameters that are not empirical. They can be deliberately
% varied and adjusted with the Catalyst Recirculation Rate (Gs) to further maximise MA Yield. 

%Sensitivity analysis is done for +/-10%, +/-20% error for the parameters b1, b2, n and n-butane:air ratio. For Inlet Temperature and 
% Pressure, absolute errors of +/-25K, +/-50K and +/-10000Pa, +/-20000Pa are used respectively. 

%For sensitivity analysis, with 5 varying values of the parameter (e.g. Original Parameter, +10% error, -10% error, +20% error, -20% error),
% the profiles for Optimum Gs (e.g. Flow of MA vs Riser Height, Temperature vs Riser Height etc.) are all plotted using the SAME Optimum Gs 
% found with the current parameters. The central idea here is that sensitivity analysis aims to see how robust the system is to errors in 
% the parameter when the reactor design is already in place. For flowrates of reagents vs Riser Height, only Maleic Anhydride and n-butane 
% are shown for sensitivity analysis. Plots for the other reagents are shown when sensitivity analysis is not selected. 

%ode23 is used for the purpose of this reactor design. 

%Integral function is used to numerically integrate the heat capacities for Adiabatic system. This takes up a huge bulk of the time 
% in which the code is run for. Manual integration (Not Recommended) can be used instead (which is currently commented out). If there 
% is a need to save time, comment out Lines 693-708 and uncomment Lines 711-715.

%% INLET PARAMETERS/ INITIAL CONDITIONS
%'o' stands for the original parameter without sensitivity. Parameters without
%'o' are those that are not required to be changed for sensitivity analyis. 
o.nBut_i = 0.01;        %Inlet n-Butane flowrate [kmol/s]
nMA_i = 0;              %Inlet Maleic anhydride flowrate [kmol/s]
nCO2_i = 0;             %Inlet Carbon Dioxide flowrate [kmol/s]
nCO_i = 0;              %Inlet Carbon Monoxide flowrate [kmol/s]
nH2O_i = 0;             %Inlet Steam flowrate [kmol/s]
o.nO2_i = 0.2;          %Inlet Oxygen flowrate [kmol/s]
o.nN2_i = 0.8;          %Inlet Nitrogen flowrate [kmol/s]
o.T_i = 720;            %Inlet Temperature [K] 
o.P_i = 410000;         %Inlet Pressure [Pa] 
psi1_i = 0;             %Inlet Extent of Raction 1 [kmol/s]. This will be 0 at the inlet as no reaction would have taken place.
c.d = 1.6;              %Riser Diameter [m]
c.zExit = 18;           %Riser Height [m] (i.e. value of z at exit) 
c.Ac = (pi.*c.d.^2)./4; %Riser Cross Sectional Area [m^2]
%%% END FOR INLET PARAMETERS/ INITIAL CONDITIONS %%%%

%% DEFINING THE CONSTANTS USED 
% 'c.' creates a structural array for all the constants, such that all the constants
%can be called altogether into the ODE function later through only one variable 'c' 
c.g = 9.81;             %Gravitational acceleration [m/s^2]
c.R = 8314.41;          %Molar Gas constant [J/kmolK]
c.pc = 1500;            %Catalyst Density [kg/m^3]
c.b1 = 429;             %Kinetic reaction parameter [m^3/kmol]
c.b2 = 2024;            %Kinetic reaction parameter [m^3/kmol]
c.n = 0.282;            %Kinetic reaction parameter 

%Pre-exponential values for each reaction
c.lnk0_Rxn1 = 2.901; c.lnk0_Rxn2 = -3.443; c.lnk0_Rxn3 = -3.938; c.lnk0_Rxn4 = 2.265; c.lnk0_Rxn5 = 3.284;

%E/R for each reaction [K]
c.ER_Rxn1 = 10672; c.ER_Rxn2 = 7012; c.ER_Rxn3 = 7012; c.ER_Rxn4 = 7012; c.ER_Rxn5 = 7012;

%Standard Enthalpy of Reaction at 298K [kJ/kmol]
c.H_Rxn1 = -1260600; c.H_Rxn2= -1526200; c.H_Rxn3= -2658600; c.H_Rxn4=-265600; c.H_Rxn5=-1398000;

%Shomate Equation Coefficients
c.A_But=0.10698; c.A_MA=-6.57385; c.A_CO2=24.99735; c.A_CO=25.56759; c.A_H2O=30.09200; c.A_O2=23.4333; c.A_N2=19.50583;
c.B_But=386.38764; c.B_MA=420.90139; c.B_CO2=55.18696; c.B_CO=6.09613; c.B_H2O=6.832514; c.B_O2=21.3322; c.B_N2=19.88705;
c.C_But=-198.80145; c.C_MA=-359.58974; c.C_CO2=-33.69137; c.C_CO=4.054656; c.C_H2O=6.793435;c.C_O2=-13.5029; c.C_N2=-8.598535;
c.D_But=39.79113; c.D_MA=116.48471; c.D_CO2=7.948387; c.D_CO=-2.671301; c.D_H2O=-2.534480; c.D_O2=3.5781; c.D_N2=1.369784;
c.E_But=0; c.E_MA=0; c.E_CO2=-0.136638; c.E_CO=0.131021; c.E_H2O=0.082139; c.E_O2=0; c.E_N2=0.527601;
%%% END FOR DEFINING CONSTANTS USED %%%

%% PROMPT TO SELECT ADIABATIC OR ISOTHERMAL SYSTEM 
global s   %"s" stands for system.
[s,v] = listdlg('PromptString','Adiabatic or Isothermal?','SelectionMode','single',...
                'ListString',{'Adiabatic', 'Isothermal'}); 
%Creates a list selection dialogue box, with integer outputs into variables "s" and "v".
%s: Output for the selection of 'Adiabatic' or 'Isothermal'. If 'Adiabatic' is 
    %selected, s=1. If 'Isothermal', s=2. 
%v: Output for 'Ok' or 'Cancel' in the dialogue box. If a system is selected 
    %and 'Ok' is selected, v=1. If 'Cancel' is selected. v=2.

if v==0     %If 'Cancel' (Invalid) Option selected.
    h=msgbox('Invalid Option. Please choose a system and run again','Error','error');
    return  %Ends the function so that the user can re-run the code with a valid option.
end

if s==1         %If 'Adiabatic' is selected.
    disp ('For Adiabatic system...')
elseif s==2     %If 'Isothermal' is selected.
    disp ('For Isothermal system...')
end
%%% END FOR ADIABATIC/ ISOTHERMAL PROMPT %%%

%% PROMPT TO CHOOSE SENSITIVITY ANALYSIS OR NOT
choice = questdlg('Would you like to do sensitivity analysis?', ...
	'Please select an option','Yes','No','Yes');
%Creates a question dialogue prompt, with character outputs 'Yes' and 'No' into 
% the character array "choice". 

%Handling response. The character output then needs to be switched to an integer 
% output so that the latter can be used in an "if" condition.
switch choice
    case 'Yes'
        sensitivity = 1;  %If 'Yes' is selected, the integer output for sensitivity = 1.
    case 'No'
        sensitivity = 2;  %If 'No' is selected, the integer output for sensitivity = 2.
end

if sensitivity == 1     %If sensitivity analysis selected.
    [p,u] = listdlg('PromptString','Which variable to change?',...
                'SelectionMode','single',...
                'ListString',{'b1', 'b2', 'n', 'Catalyst Density', 'n-butane:air ratio',...
                'Inlet Temperature','Inlet Pressure'});
%Creates a list selection dialogue box, with integer outputs into variables "p" and "u".
%"p" stands for parameter.
% p: Output for the selection of the parameter to do sensitivity analysis on. 
%   The integer outputs produced are based on the order of the 7 variables in ListString 
%   (e.g. If 'b1' is selected, p=1. If 'n' is selected, p=3. If 'Inlet Pressure', p=7 etc.)
% u: Output for 'Ok' or 'Cancel' in the dialogue box. If a system is selected 
%   and 'Ok' is selected, u=1. If 'Cancel' is selected. u=2. 
     disp('Please wait. Calculation in progress...')
     if u==0    %If 'Cancel' Option selected, NO sensitivtity analysis.
         p=8;   %p will be = 8 if sensitvity is not selected.
         disp('No sensitivity selected.')
     end
elseif sensitivity == 2   %If NO sensitivity analysis selected.
    p=8; 
    disp('Sensitivity not selected.')
    disp('Please wait. Calculation in progress...')
end
%%% END FOR SENSITIVITY ANALYSIS PROMPT %%%

%% PROMPT FOR FURTHER SIMULATION (ONLY IF SENSITIVITY OF 'n-butane:air ratio','Inlet Temperature','Inlet Pressure' SELECTED)
if p==5 || p==6 || p==7  %Integer outputs into "p" if 'n-But:air ratio', 'Inlet T' or 'Inlet P' selected. 
choice2 = questdlg('Would you like to do Further Simulation/ Optimisation?',...
	'Please select an option', 'Yes','No','Yes');
%Similar to the question dialogue box used in sensitivity analysis prompt.
switch choice2
    case 'Yes'
        FurtherSimulation = 1;
        disp('Further Optimisation plots will be displayed.')
    case 'No'
        FurtherSimulation = 2;
        disp('Further Optimisation not selected.')
end
%Here, the character output in character array "choice" is switched to an
% integer output in "FurtherSimulation".
end
%%% END FOR FURTHER SIMULATION PROMPT %%%

%% DEFINING VARIABLES IN THE LOOPS AS MATRICES 
MA_o=[];    %Outlet MA Flowrate [kmol/s]
T=[];       %Outlet Temperature [K]
B_c=[];     %n-Butane Conversion [%]
MA_s=[];    %MA Selectivity [%]
%These matrices are used in the loops and have to be defined before the outer 'N' loop below
% as writing them inside the loops will cause the matrix to re-define itself into an empty 
% matrix each time a new iteration is executed. Thus, in order to prevent the data in the 
% matrices from overwriting at the start of each 'N' iteration, they are first defined here.  
%%% END FOR DEFINING VARIABLES AS MATRICES SECTION %%

%% SETTING UP OUTER 'N' LOOP AND VARYING SELECTED PARAMETERS FOR SENSITIVITY ANALYSIS 
varypercent = [1 1.1 0.9 1.2 0.8]; %Creates a matrix of +/- 10%, +/- 20% error for the Parameter
varyabsoluteT = [0 25 -25 50 -50]; %Creates a matrix of +/- 25K, +/- 50K absolute error for Inlet Temperature
varyabsoluteP = [0 10000 -10000 20000 -20000]; %Creates a matrix of +/- 10000Pa, +/- 20000Pa absolute error for Inlet Pressure

%Inlet Temperature and Pressure will be varied by absolute errors in the sensitivity analysis. 
% Thus, these absolute errors will be added to the orginal inlet Temperature. 

%For the other parameters, they will be varied by percentage errors. Thus, these values in 
% the matrix will be multiplied to the original parameters. (e.g. When +10% error of b1 is 
% desired, b1(new) = b1(original)*1.1)

%Note that the first values in the matrices above are 1 or 0. Thus, the % or absolute error in the 
% first iteration of the loop will be *1 or +0 respectively, i.e. when the parameter is original.

b1=429; b2=2024; n=0.282; pc = 1500;    
nBut_i=o.nBut_i; nO2_i=o.nO2_i; nN2_i=o.nN2_i; T_i=o.T_i; P_i=o.P_i;
%Constant parameters re-defined before loop (without 'c.' in front). Explained below.
%Inlet conditions re-defined (without 'o.' in front). Explained below.

for N=1:5   %OUTER LOOP - for 5 different values of the parameter chosen for sensitivity analyis.
    if p==1     %Vary b1
    c.b1 = varypercent(N)*b1;   %The parameters called into the ODE Function are those with 
    % 'c.' in front. When changing the parameter for sensitivity analysis, the error is 
    % multiplied by the parameters without 'c.' (which remain constant) in Line 175. This 
    % ensures the adjustments with the errors are not cumulative from the previous iteration. 
    elseif p==2     %Vary b2
    c.b2 = varypercent(N)*b2;    
    elseif p==3     %Vary n
     c.n = varypercent(N)*n;  
    elseif p==4     %Vary pc
    c.pc = varypercent(N)*pc;   
    elseif p==5     %Vary n-butane:air ratio
    %For inlet conditions however, the inlet conditions called into the ODE function are
    % those WITHOUT 'o.'. Thus, the error is multipled by/ added to the parameter with
    % 'o.' (which remain constant). Again, this ensures the adjustments with the errors are
    % not cumulative from the previous iteration. 
    
    %Total flowrate n_t needs to be kept constant (1.01 kmol/s in this case). Thus, the 
    % n-butane:air ratio cannot be varied by solely varying the n-butane flowrate. Flowrates 
    % of air (Nitrogen and Oxygen) also needs to be decreased by a certain amount to keep n_t 
    % constant. The mathematic expressions for the flowrates to attain the desired percentage 
    % changes of n-butane:air ratios are worked out and presented as follows.
    nBut_i = (1.01*varypercent(N)*o.nBut_i)./(1.01+(o.nBut_i)*(varypercent(N)-1));
    nO2_i = o.nO2_i + (o.nBut_i-nBut_i)*0.2;
    nN2_i = o.nN2_i + (o.nBut_i-nBut_i)*0.8;
    %Flowrate of air needs to be decreased/ increased by the amount in which nBut_i is 
    % increased/ decreased. To maintain the proportions in air, nO2_i is increased/decreased
    % by 20% of this amount while nN2_i by 80%.
    elseif p==6     %Vary inlet Temperature
     T_i = (varyabsoluteT(N)+o.T_i);
    elseif p==7     %Vary inlet Pressure 
     P_i = (varyabsoluteP(N)+o.P_i); 
    else    %No sensitivity. All parameters as original.
    end
    %%% END FOR VARYING SELECTED PARAMETERS FOR SENSITIVITY ANALYSIS %%%
    
%% SETTING UP INNER 'Gs' LOOP TO RUN ODE WITH VARYING VALUES OF CATALYST RECIRCULATION RATE GS
z = 0:0.01:c.zExit; %For each value of Gs, the dependent variables are calculated at z intervals of 0.01

%Adjusting the range of Gs used in the loop for the 2 different systems.
if s==1     %If Adiabatic system
    Gs_low = 1; Gs_step = 10; Gs_final = 1001;
elseif s==2     %If Isothermal system
    Gs_low=1;Gs_step = 50; Gs_final = 5001;
end
%For Adiabatic, max MA outlet flowrate is achieved at a lower optimum Gs, hence the lower Gs_final. 
%To increase the efficiency and speed of the code, the Gs_step is such that only 101 iterations
% are executed in the inner Gs loop. The values obtained from these 101 iterations of Gs are then
% used to plot curves that follows. The curves will not be the smoothest, but smooth enough to 
% observe the trends. The exact value (to the nearest 1) of the Optimum Gs will be found later.
for Gs = Gs_low:Gs_step:Gs_final
 DV_i = [nBut_i nMA_i nCO2_i nCO2_i nCO_i nO2_i nN2_i T_i P_i psi1_i]; %Initial Parameters of Dependent Variables (DVs)
 [zSol,DVsol] = ode23(@(z,DV) ODESolver(z,DV,Gs,c), z,DV_i);  %Solve by the ode23 integrator.
        %Two input arguments, z and DV. Pass in extra parameters by defining them outside the
        % ODESolver Function and passing them in when specifying the function handle.
        %A resulting 1801x10 matrix of "DVsol" is produced: 1801 rows (1801 steps of z) and 10 
        % columns (10 DVs). Thus 1st row represents the inlet value of the DV and 
        % 1801th (or length(z))th row represents the outlet value of the DV.
 MA_outlet.N = DVsol(length(z),2);  %Outlet Flowrates of the MA product for each iteration of Gs. 
 T_outlet.N = DVsol(length(z),8);   %Outlet Temperatures for each iteration of Gs.
 But_conv.N = ((DVsol(1,1)-DVsol(length(z),1))./(DVsol(1,1))).*100; %n-Butane Conversion for each iteration of Gs. 
 MA_selectivity.N = ((DVsol(length(z),2))./(DVsol(1,1)-DVsol(length(z),1))).*100; %MA Selectivity for each iteration of Gs.  
 
 nn=((Gs+Gs_step-1)./Gs_step);  %nn gives what iteration the Gs (1st, 2nd, 3rd... etc.) is on. Used for indexing later.
 
 Gsrange(nn) = Gs;  %Stores the values of Gs iterations into a single row matrix.
 
 MA_o(nn,N) = MA_outlet.N;  %For Outlet MA Flowrate.  
 T(nn,N) = T_outlet.N;  %For Outlet Temperature. 
 B_c(nn,N) = But_conv.N;    %For n-Butane Conversion. 
 MA_s(nn,N) = MA_selectivity.N; %For MA Selectivity.
 %The 4 matrices above are created into a 101x5 matrix, which stores the respective data necessary.
 % Using indexing, the 'nn'th row represent what iteration Gs is on, and 'N'th column what iteration 
 % N is on. For e.g., for sensitivity of b1 in adiabatic system, MA_o(5,2) represents the outlet MA 
 % flowrate at the 5th iteration of Gs (When Gs=55) and N=2 (When c.b1= b1*varypercent(2)= 429*1.1).
end
%%% END FOR INNER GS LOOP %%%

%% PLOTTING REQURIED GRAPHS FROM THE MATRICES GENERATED IN THE INNER GS LOOP
%Overall Performance as a function of Catalyst Recirculation Rate Gs. These
% 'plot' functions are still within the outer N loop. Hence, "hold on" is
% used before each plot to ensure the curve for each N iteration (each varied
% parameter in question) stays in the same figure for easy comparison. 
%As shown below, the plots are given variables (e.g. FigOne, FigTwo etc.).
% These are useful for the legends which will be shown later. 
figure (1) 
x = Gsrange; y = MA_o(:,N); %x values are the values of Gs iterations and y values are 
    %the resulting MA outlet flowrates achieved (at a certain N).
idxmax = find(y == max(y));     %Finds the index of the maximum Outlet MA Flowrate. 
hold on
FigOne(N)= plot(x,y, '-h','MarkerIndices',[idxmax],...
    'MarkerEdgeColor','red', 'MarkerFaceColor', 'red',...
    'MarkerSize',5); %Plots the graph of Outlet MA flowrate vs Gs and displays a red marker at the
    % maximum point. After one loop of N, it holds on and plots for the next column of MA outlet flowrates. 
xlabel('Catalyst Recirculation Rate Gs [kgm^{2}s^{-1}]')
ylabel('Outlet Flowrate of MA [kmols^{-1}]')
title('Outlet Flowrate of MA vs. Catalyst Recirculation Rate')
xlim([0 Gs_final])

figure(2)
hold on
FigTwo(N)= plot(Gsrange, T(:,N));   %Plots the graph of Outlet Temperature vs Gs.
xlabel('Catalyst Recirculation Rate Gs [kgm^{2}s^{-1}]')
ylabel('Outlet Temperature [K]')
title('Outlet Temperature vs. Catalyst Recirculation Rate')
xlim([0 Gs_final])

figure (3) 
hold on
FigThree(N)= plot(Gsrange,B_c(:,N));    %Plots the graph of n-Butane Conversion vs Gs.
xlabel('Catalyst Recirculation Rate Gs [kgm^{2}s^{-1}]')
ylabel('n-Butane Conversion [%]')
xlim([0 Gs_final])
ylim([0 100])
title('n-Butane Conversion vs. Catalyst Recirculation Rate')

figure (4) 
hold on
FigFour(N)= plot(Gsrange,MA_s(:,N));    %Plots the graph of MA Selectivity vs Gs.
xlabel('Catalyst Recirculation Rate Gs [kgm^{2}s^{-1}]')
ylabel('MA Selectivity [%]')
xlim([0 Gs_final])
ylim([0 100])
title('MA Selectivity vs. Catalyst Recirculation Rate')

figure (5)
hold on
FigFive(N)= plot(B_c(:,N),MA_s(:,N));   %Plots the graph of MA Selectivity vs n-Butane Conversion
xlabel('n-Butane conversion [%]')
ylabel('MA Selectivity [%]')
ylim([0 100])
title('MA Selectivity vs. n-Butane Conversion')
%%% END FOR GRAPHS IN INNER GS LOOP/ OVERALL PERFORMANCE GRAPHS %%%

%% 	FINDING THE EXACT OPTIMUM GS TO BE USED IN REACTOR (WHEN MAXIMUM OUTLET MA FLOWRATE IS ACHIEVED)
[maxMA, location] = max(MA_o(:,N)); %Gives the max outlet MA flowrate in "maxMA" and the 
        % row of matrix at which this is achieved in "location" for each N iteration.
Row_maxMA(N) = location;  %Stores the "location" for each N iteration.

%Previously, a step of 10 was used for the Gs loop. Thus, the optimum Gs 
% that will be found will only be to the nearest 10 (e.g. 241, 251, 261 etc.). 
% In order to find the precise value of the optimum Gs, the ODEFunction is ran 
% again in steps of 1 within a shorter range where the optimum Gs is definitely 
% present. This range will be: (Optimum Gs found to the nearest 10) +/- Gs_step.
%That way, the code is ran more efficiently -- less time taken but with
% precision still present. 

MA_oexact=[];
    OptimumGs_Rounded(N) = (Row_maxMA(N)-1)*Gs_step + 1; %Converting the numbered row into the value 
                % of Optimum Gs (to the nearest 10) [kgcat/m^2*s].
    for Gs = (OptimumGs_Rounded(N)-Gs_step):1:(OptimumGs_Rounded(N)+Gs_step) %Running the ODEFunction 
                % again in steps of 1 within a shorter range for a precise value of Optimum Gs.
DV_i = [nBut_i nMA_i nCO2_i nCO2_i nCO_i nO2_i nN2_i T_i P_i psi1_i]; %Initial parameters of dependent variables
[zSolb,DVsolb] = ode23(@(z,DV) ODESolver(z,DV,Gs,c), z,DV_i);
MA_outletexact = DVsolb(length(z),2);  %Outlet Flowrates of the MA product for each iteration of Gs
MA_oexact = [MA_oexact;MA_outletexact,Gs];  %Creates a 2 column matrix to store: Outlet MA Flowrate and Gs at which this is achieved.
[maxMAexact,location] = max(MA_oexact(:,1));  %Gives the max outlet MA flowrate in "maxMAexact" and the 
                % row of matrix at which this is achieved in "location". 
maxMA_Exact(N) = maxMAexact;   %Stores the max MA outlet flowrate for each iteration of N
OptimumGs_Exact(N) = MA_oexact(location,2);  %The precise value of the Optimum Gs can then be found. Indexing
                % is used to store the exact Optimum Gs for each iteration of N.
    end
OptimumRecirculationRateGs = OptimumGs_Exact(1); %Selects the first optimum Gs (i.e. found with original parameters)
%%% END FOR FINDING OPTIMUM GS %%%

%% PROFILES FOR OPTIMUM CATALYST RECIRCULATION RATE GS
%The following profiles are plotted using the first Optimum Gs value that was
% found in the previous section (i.e. with original parameters), even for
% the remaining 4 iterations of N. The central idea here is that sensitivity
% analysis aims to see how robust the system is to errors when already in place. 
 
%Now, the aim is to use a fixed Gs to attain the riser profiles with height.
% Thus, an outer loop of the Optimum Gs is created with an inner loop executed
% over the range of z (riser height). 

for Gs=OptimumRecirculationRateGs       %As explained above.
DV_i = [nBut_i nMA_i nCO2_i nCO2_i nCO_i nO2_i nN2_i T_i P_i psi1_i]; %Initial Parameters
 [zSolc,DVsolc] = ode23(@(z,DV) ODESolver(z,DV,Gs,c), z,DV_i);
DV_3D(:,:,N) = DVsolc; %Stores the data for the Dependent Variables generated here 
            % in a 3-D array. 5 iterations of 'N' produces 5 pages, where each page 
            % gives the 1801x10 double matrix of DVs. These data are needed later. 
            
One_voidage=[];     %Defines 1-voidage as a matrix.
CO_s =[];           %Defines the CO selectivity as a matrix.
CO_over_CO2 = [];   %Defines CO/CO2 as a matrix.
%These 3 matrices will fill up below, and the values will be used to plot the
% curves for each 'N' iteration. Since they are defined in the 'N' loop, they will be 
% re-defined as empty matrices and be overwritten at the start of each 'N' iteration. 
% However, the intended purpose here is only to plot the curves for comparison, it is 
% not necessary to keep the data and this is thus neglected.
 for i=1:length(z)      %Creating an inner loop for the range of z. 
     n_t = (DVsolc(i,1)+DVsolc(i,2)+DVsolc(i,3)+DVsolc(i,4)+DVsolc(i,5)+DVsolc(i,6)+DVsolc(i,7));   %Total Flowrate [kmol/s] at a specific 'z'.
     v_t = (n_t*c.R.*DVsolc(i,8))./DVsolc(i,9);     %Volumetric Flowrate [m^3/s] at a specific 'z'.
     U_0 = v_t./c.Ac;   %Riser superficial gas velocity [m/s] at a specific 'z'.
     One_voidage(i) = 1-(U_0*c.pc)./(2*Gs+(U_0*c.pc));  %Matrix of 1-voidage for each iteration of 'z'.
     CO_s(i) = ((DVsolc(i,4))/(4*(DVsolc(1,1)-DVsolc(i,1))))*100;   %Matrix of CO selectivity [%] for each iteration of 'z'.
     CO_over_CO2(i) = (DVsolc(i,4))./(DVsolc(i,3)); %Matrix of CO/CO2 flowrate for each iteration of 'z'
 end
end

%Like before, "hold on" is used  to ensure the curve for each N iteration (each 
% varied parameter in question) stays in the same figure for easy comparison. 
%% MA & n-Butane Flowrates vs. Riser Height (For comparison in Sensitivity Analysis)
figure(6)
FigSix(N)=plot(zSolc,DVsolc(:,2));  %Plots graph of MA Flowrate vs Riser height
hold on
title('Flowrate of MA vs. Riser height')
xlabel('Riser Height [m]')
ylabel('Flowrate of MA [kmol/s]')

figure(7)
FigSeven(N)=plot(zSolc,DVsolc(:,1));    %Plots graph of n-Butane Flowrate vs Riser height
hold on
title('Flowrate of n-Butane vs. Riser height')
xlabel('Riser Height [m]')
ylabel('Flowrate of n-Butane [kmol/s]')

%% Temperature and Pressure Profiles
figure(8) 
FigEight(N)=plot(zSolc,DVsolc(:,8));    %Plots graph of Temperature vs Riser height
hold on
title('Temperature vs. Riser height')
xlabel('Riser Height [m]')
ylabel('Temperature [K]')

figure(9)
FigNine(N)=plot(zSolc,DVsolc(:,9));     %Plots graph of Pressure vs Riser height
hold on
title('Pressure vs. Riser height')
xlabel('Riser Height [m]')
ylabel('Pressure [Pa]')

%% 1-voidage vs. Riser Height
figure(10)
FigTen(N)=plot(zSolc,One_voidage);  %Plots graph of 1-voidage vs Riser height
hold on
title('1-\epsilon_a_v_e vs. Riser height')
xlabel('Riser Height [m]')
ylabel('1-\epsilon_a_v_e')

%% CO Selectivity and CO/CO2 vs. Extent of Reaction 1
figure(11)
FigEleven(N)=plot(DVsolc(:,10),CO_s); %Plots graph of CO Selectivity vs Extent of Reaction 1
hold on
title('CO Selectivity vs. Extent of Reaction 1')
xlabel('Extent of Reaction 1[kmol/s]')
ylabel('CO Selectivity [%]')

figure(12)
FigTwelve(N)=plot(DVsolc(:,10),CO_over_CO2); %Plots graph of CO/CO2 vs Extent of Reaction 1
hold on
title('CO/CO2 vs. Extent of Reaction 1')
xlabel('Extent of Reaction 1[kmol/s]')
ylabel('CO/CO2 Flowrate')

%% BREAKING THE LOOP IF NO SENSITIVITY ANALYSIS SELECTED
if p==8   
    format shortE
    MaxMAOutletFlowrate = maxMA_Exact   %Displays the max MA Outlet Flowrate from the first 'N' iteration.
    disp('[kmol/s]')        
figure(13)  %Graph of the flowrates of all reagents combined into one figure for easier viewing. 
plot(zSolc,DVsolc(:,2))           %Plots graph of MA Flowrate vs Riser height
hold on, plot(zSolc,DVsolc(:,1))  %Plots graph of n-Butane Flowrate vs Riser height
hold on, plot(zSolc,DVsolc(:,3))  %Plots graph of CO2 Flowrate vs Riser height
hold on, plot(zSolc,DVsolc(:,5))  %Plots graph of H2O Flowrate vs Riser height
hold on, plot(zSolc,DVsolc(:,4))  %Plots graph of CO Flowrate vs Riser height
title('Flowrates vs. Riser height')
xlabel('Riser Height [m]')
ylabel('Flowrate [kmol/s]')
legend ('MA','n-Butane','CO2','H2O','CO','Location','northwest')
    break   %Breaks the loop after the first 'N' iteration (when percentage error is 1 and 
            %absolute error is 0, i.e. all parameters original).
end
end    %End for outermost 'N' loop% 
OptimumRecirculationRateGs  %Displays in the Command Window the Optimum Gs to be used for the reactor design.
disp('[kgm^2/ s]')
%%% END FOR BREAKING LOOP IF NO SENSITIVITY ANALYSIS %%%
%%% END FOR PROFILES FOR OPTIMUM GS %%%

%% QUANTIFYING PERCENTAGE CHANGE OF MA YIELD WITH %ERROR/ ABSOLUTE ERROR IN PARAMETER
if p==8
    return  %Returns control as the following is not necessary if no sensitivity analysis selected.
end

%MA Yield defined as: (Outlet Flowrate of MA)/(Inlet Flowrate of n-Butane)%
%First, finding the Outlet MA Flowrate for 5 different values of the parameter 
% chosen for sensitivity analyis from the 3-D array "DV_3D" created in Line 359.
MAoutlet3D = [];
MAoutlet3D = [MAoutlet3D; DV_3D(1801,2,:)]; %1x1x5 Matrix of the outlet MA Flowrates for the 5 parameters 
MAoutlet_sensitivity = squeeze(MAoutlet3D); %Removes singleton dimensions to make it 2D

%Second, finding the Inlet n-But Flowrate for 5 different values of the parameter.
Butinlet3D = [];
Butinlet3D = [Butinlet3D; DV_3D(1,1,:)]; %1x1x5 Matrix of the inlet But Flowrates for the 5 parameters
Butinlet_sensitivity = squeeze(Butinlet3D); %Removes singleton dimensions to make it 2D

for N=1:5
%Next, finding the MA Yield for each of the 5 different values of the parameter
    MAYield(N) = (MAoutlet_sensitivity(N))./(Butinlet_sensitivity(N));
    
%Finally, finding the percentage change in MA Yield
    MAYieldChange(N) = ((MAYield(N)-MAYield(1))./(MAYield(1)))*100;
    percentageerror(N) = [(varypercent(N)-1)*100];  %Converting [1 1.1 0.9 1.2 0.8] into %error (i.e. 0%,+10%,-10%... etc.) 
end

if p<=5    %If first 5 parameters selected for sensitivity. Error for these parameters are %error.
    xerror = percentageerror;
elseif p==6    %If T selected. Error for T is in absolute.
    xerror = varyabsoluteT;
elseif p==7    %If P selected. Error in absolute.var
    xerror = varyabsoluteP;
end
fit_xerror = linspace(min(xerror),max(xerror),100);
fit_MAYieldChange = interp1(xerror,MAYieldChange,fit_xerror,'spline');
%As values of the errors are not in ascending or descending order (e.g. 0,10, -10, 20, -20), 
% the plots will be zigzagged as Matlab joins the first value of matrix to the second point, 
% and second to the third etc. Thus, linspace generates 100 evenly spaced points between the 
% lowest to the highest value of the error. interp1 returns interpolated values of 
% MAYieldChange between each data point to generate a smooth curve.  

figure(13)
scatter(xerror,MAYieldChange,25,'filled')   %Plots the 5 data points used to generate the curve.
hold on
FigThirteen = plot(fit_xerror,fit_MAYieldChange);   %Plots graph of %Change in MA Yield vs Error in Parameter.
grid on
ax = gca;   %current axes
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';    %Controls axis location 
xlabeltext = {'Error in Kinetic Parameter b1 [%]', 'Error in Kinetic Parameter b2 [%]',...
    'Error in Kinetic Parameter n [%]', 'Error in Catalyst Density [%]',...
    'Error in n-But:air ratio [%]','Error in Temperature [K]', 'Error in Pressure [Pa]'};
xlabel(xlabeltext{p})   %Changes the x-axis label depending on the parameter selected for sensitivity.
ylabel('Change in MA Yield [%]')
title('% Change in MA Yield vs. Error in Parameter (with current design)')
title('Extent in which MA Yield will be affected by the Errors in Parameter')
%%% END FOR PERCENTAGE CHANGE OF MA YIELD VS PERCENTAGE/ABSOLUTE ERROR IN PARAMETER %%%

%% LEGENDS FOR PLOTS 
%Since different parameters can be selected for sensitivity analysis, the
% code needs to be written in such a way that the legend changes based on
% the parameters selected for sensitivity analysis.
legendtext1 = {'Kinetic Parameter b1', 'Kinetic Parameter b2','Kinetic Parameter n','Catalyst density',...
    'n-But:air ratio','Inlet Temperature','Inlet Pressure'};    %Order of parameters in the cell array 
                                %similar to order of the output for p.
legendtext2 = {'with 10% error', 'with -10% error', 'with 20% error', 'with -20% error',...
    'of +25K', 'of -25K', 'of +50K', 'of -50K',...
    'of +10kPa', 'of -10kPa','of +20kPa','of -20kPa'};   %3 sets of 4 character vectors - 1st set for 
                                % varypercent, 2nd set for varyabsoluteT, 3rd set for varyabsoluteP.

if p<=5
    L = 1:4;    %Ensures 1st set of legendtext2 is used.
elseif p==6
    L = 5:8;    %Ensures 2nd set of legendtext2 is used.
elseif p==7
    L = 9:12;   %Ensures 3rd set of legendtext2 is used.
end

%Combining the two character vectors "legendtext1" and "legendtext2", which both changes 
% based on the parameter selected.
firstlegend = 'Original Parameters';
secondlegend = sprintf('%s %s',legendtext1{p},legendtext2{L(1)});  
thirdlegend = sprintf('%s %s',legendtext1{p},legendtext2{L(2)});
fourthlegend = sprintf('%s %s',legendtext1{p},legendtext2{L(3)});
fitfhlegend = sprintf('%s %s',legendtext1{p},legendtext2{L(4)});

%Now, putting the 5 legends into each Figure.
if s==1
legend(FigOne, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend)
elseif s==2
legend(FigOne, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend,'Location','southeast')
end
legend(FigTwo, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend,'Location','southeast')
legend(FigThree, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend,'Location','southeast')
legend(FigFour, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend)
legend(FigFive, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend)
legend(FigSix, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend,'Location','northwest')
legend(FigSeven, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend)
legend(FigEight, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend,'Location','southeast')
legend(FigNine, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend)
legend(FigTen, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend)
legend(FigEleven, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend)
legend(FigTwelve, firstlegend, secondlegend, thirdlegend, fourthlegend, fitfhlegend)
%%% END FOR LEGENDS %%%

%% FURTHER SIMULATION: FURTHER OPTIMISATION OF MA YIELD
%The next step is to see how to further optimise MA Yield. Parameters like Initial Feed Ratio, 
% Inlet T and Inlet P can be adjusted to give a higher yield. However, these max MA Yields 
% are attained at different Optimum Gs and thus this section aims to quantify the max MA Yield 
% possible over a range of the 3 parameters mentioned above. Also, the optimum Gs at which
% these maximum yields are achieved. 
for N=1:5
    MaxMAYield_Optimisation(N) = (maxMA_Exact(N)/Butinlet_sensitivity(N))*100;
    %Quantifying Maximum MA Yield (when optimum Gs is used) = (Max MA
    % Outlet)/(Inlet n-Butane) * 100%
end
ButAirRatio_OptimisationRange = (o.nBut_i/(o.nO2_i+o.nN2_i))*(varypercent);
T_OptimisationRange = o.T_i+varyabsoluteT;
P_OptimisationRange = o.P_i+varyabsoluteP;
%Above vectors are the 5 values of each parameter for which the Optimum Gs and max MA Yield 
% will be determined for.
%% Optimising Initial Feed Ratio/ Inlet Temperature/ Inlet Pressure
if p==5
    xRangeOptimisation = ButAirRatio_OptimisationRange;
elseif p==6
    xRangeOptimisation = T_OptimisationRange;
elseif p==7
    xRangeOptimisation = P_OptimisationRange;
else
    return  %If sensitivity of other parameters not selected, code stops running.
end
if FurtherSimulation==2
    return   %If Further Simulation is not selected, then code also stops running.
end

fit_xRangeOptimisation = linspace(min(xRangeOptimisation),max(xRangeOptimisation),100);
fit_OptimumGs_Optimisation = interp1(xRangeOptimisation,OptimumGs_Exact,fit_xRangeOptimisation,'spline');
fit_MaxMAYield_Optimisation = interp1(xRangeOptimisation,MaxMAYield_Optimisation,fit_xRangeOptimisation,'spline');
%Again, since only 5 points are used, the curve will not be smooth.
% linspace generates 100 evenly spaced points between the lowest to the highest value of the parameter. 
% interp1 is used to return interpolated y-values between each data point to generate a smooth curve.  

figure(14)
%Graph for Optimum Gs vs. Range of Parameter
subplot(1,2,1)
plot(fit_xRangeOptimisation,fit_OptimumGs_Optimisation)
hold on
CurrentParameter = plot(xRangeOptimisation(1),OptimumGs_Exact(1),'r.','MarkerSize',15); %Displays the point
                        % which current design parameters are on.
legend([CurrentParameter], {'Current Parameter'},'Location','northeast')
%For axis label (x-axis label changes based on parameter chosen):
xlabeltext2 = {'','','','','n-But:air ratio','Inlet Temperature [K]', 'Inlet Pressure [Pa]'};
xlabel(xlabeltext2{p})
ylabel('Optimum Gs [kgm^{2}s^{-1}]') 
%For title (which changes based on parameter chosen):
titletext = {'','','','','Optimum Catalyst Recirculation Rate Gs vs. n-But:air ratio',...
    'Optimum Catalyst Recirculation Rate Gs vs. Inlet Temperature',...
    'Optimum Catalyst Recirculation Rate Gs vs. Inlet Pressure'};
title(titletext{p})

%Graph for Max MA Yield Possible vs. Range of Parameter
subplot(1,2,2)
plot(fit_xRangeOptimisation, fit_MaxMAYield_Optimisation)
hold on
CurrentParameter = plot(xRangeOptimisation(1),MaxMAYield_Optimisation(1),'r.','MarkerSize',15); %Displays the point
                        % which current design parameters are on.
legend([CurrentParameter], {'Current Parameter'},'Location','southeast')
xlabel(xlabeltext2{p})
ylabel('Max MA Yield [%]')
titletext2 = {'','','','','Maximum MA Yield Possible vs. n-But:air ratio',...
    'Maximum MA Yield Possible vs. Inlet Temperature', 'Maximum MA Yield Possible vs. Inlet Pressure'};
title(titletext2{p})
%%% END FOR FURTHER SIMULATION/ OPTIMISATION %%%

%% ODE SOLVER FUNCTION %%
function dDVdz = ODESolver(z,DV,Gs,c)

%% Dependent Variables - DV
nBut = DV(1);   %n-Butane flowrate [kmol/s]
nMA = DV(2);    %Maleic Anhydride flowrate [kmol/s]
nCO2 = DV(3);   %Carbon Dioxide flowrate [kmol/s]
nCO = DV(4);    %Carbon Monoxide flowrate [kmol/s]
nH2O = DV(5);   %Steam flowrate [kmol/s]
nO2 = DV(6);    %Oxygen flowrate [kmol/s]
nN2 = DV(7);    %Nitrogen flowrate [kmol/s]
T = DV(8);      %Riser Temperature [K]
P = DV(9) ;     %Riser Pressure [Pa]
psi1 = DV(10);  %Extent of Reaction 1 [kmol/s]

%% Total Flowrates in Riser
n_t = nBut+nMA+nCO2+nCO+nH2O+nO2+nN2;   %Total Molar Flowrate (kmol/s)
v_t = (n_t.*c.R.*T)./P;     %Total Volumetric Flowrate (m^3/s)

%% Riser Superficial Gas Velocity [m^2/s]
U_0 = v_t./c.Ac;

%% Axial voidage/ Riser voidage
ep = (U_0.*c.pc)./(2.*Gs + U_0.*c.pc);

%% K Constant to multiply to ODEs [kgcat/m]
K = c.Ac.*c.pc*(1-ep);

%% Concentrations of n-Butane, MA, Oxygen [kmol/m^3]
CBut = nBut./v_t;
CMal = nMA./v_t;
COxy = nO2./v_t;
 
%% Reaction constant for each reaction
k.Rxn1 = exp(c.lnk0_Rxn1 - (c.ER_Rxn1./T));
k.Rxn2 = exp(c.lnk0_Rxn2 - (c.ER_Rxn2./T)); 
k.Rxn3 = exp(c.lnk0_Rxn3 - (c.ER_Rxn3./T));
k.Rxn4 = exp(c.lnk0_Rxn4 - (c.ER_Rxn4./T));
k.Rxn5 = exp(c.lnk0_Rxn5 - (c.ER_Rxn5./T));

%% Rate Expression for each reaction [kmol/ s*kgcat]
r.Rxn1=(k.Rxn1.*(COxy.^c.n).*c.b1.*CBut)./(1+c.b1.*CBut); 
r.Rxn2=(k.Rxn2.*(COxy.^c.n).*c.b2.*CBut)./(1+c.b2.*CBut); 
r.Rxn3=(k.Rxn3.*(COxy.^c.n).*c.b2.*CBut)./(1+c.b2.*CBut); 
r.Rxn4=(k.Rxn4.*(COxy.^c.n).*CMal)./(1+c.b1.*CBut); 
r.Rxn5=(k.Rxn5.*(COxy.^c.n).*CMal)./(1+c.b1.*CBut);

%% Gas Phase Heat Capacity Cp Functions [J/ mol*K]
t = T./1000;
Cp_But = @(t) c.A_But + c.B_But.*t + c.C_But.*t.^2 + c.D_But.*t.^3 + c.E_But./(t.^2);
Cp_MA = @(t) c.A_But + c.B_But.*t + c.C_But.*t.^2 + c.D_But.*t.^3 + c.E_But./(t.^2);
Cp_CO2 = @(t) c.A_CO2 + c.B_CO2.*t + c.C_CO2.*t.^2 + c.D_CO2.*t.^3 + c.E_CO2./(t.^2);
Cp_CO = @(t) c.A_CO + c.B_CO.*t + c.C_CO.*t.^2 + c.D_CO.*t.^3 + c.E_CO./(t.^2);
Cp_H2O = @(t) c.A_H2O + c.B_H2O.*t + c.C_H2O.*t.^2 + c.D_O2.*t.^3 + c.E_H2O./(t.^2);
Cp_O2 = @(t) c.A_O2 + c.B_O2.*t + c.C_O2.*t.^2 + c.D_O2.*t.^3 + c.E_O2./(t.^2);
Cp_N2 = @(t) c.A_N2 + c.B_N2.*t + c.C_N2.*t.^2 + c.D_N2.*t.^3 + c.E_N2./(t.^2);

%% Integrals of Cp Functions (Using MATLAB Numerical Integration)
global s 
% if s==1        %Finding the Heat of Reactions is only necessary for Adiabatic system. Thus, running the
%                % integral functions for Isothermal system is unneccessary and can be ommitted to save time. 
% But_integral = integral(Cp_But, 0.298, t);      %Integrals are from 0.298K to t (Temperature of Riser/1000).
% MA_integral = integral(Cp_MA, 0.298, t);
% CO2_integral = integral(Cp_CO2, 0.298, t);
% CO_integral = integral(Cp_CO, 0.298, t);
% H2O_integral = integral(Cp_H2O, 0.298, t);
% O2_integral = integral(Cp_O2, 0.298, t);
% 
% %% Heat of Reaction for each reaction from 298K to T 
% % dH_Rxn1 = c.H_Rxn1 - But_integral - 3.5*O2_integral + MA_integral + 4*H2O_integral;
% % dH_Rxn2 = c.H_Rxn2 - But_integral - 4.5*O2_integral + 4*CO_integral + 5*H2O_integral;
% % dH_Rxn3 = c.H_Rxn3 - But_integral - 6.5*O2_integral + 4*CO2_integral + 5*H2O_integral;
% % dH_Rxn4 = c.H_Rxn4 - MA_integral - O2_integral + 4*CO_integral+ H2O_integral;
% % dH_Rxn5 = c.H_Rxn5 - MA_integral - 3*O2_integral +4*CO2_integral+ H2O_integral;
% end

%% Heat of Reaction for each reaction from 298K to T (Manual Integration) - NOT RECCOMENDED
dH_Rxn1 = c.H_Rxn1 + ((t-0.298)*(-c.A_But-3.5.*c.A_O2+c.A_MA+4.*c.A_H2O) + ((t.^2-0.298.^2)./2)*(-c.B_But-3.5.*c.B_O2+c.B_MA+4.*c.B_H2O) + ((t.^3-0.298.^3)./3)*(-c.C_But-3.5.*c.C_O2+c.C_MA+4.*c.C_H2O) + ((t.^4-0.298.^4)./4)*(-c.D_But-3.5.*c.D_O2+c.D_MA+4.*c.D_H2O) + ((t.^5-0.298.^5)./5)*(-c.E_But-3.5.*c.E_O2+c.E_MA+4.*c.E_H2O));
dH_Rxn2 = c.H_Rxn2 + ((t-0.298)*(-c.A_But-4.5.*c.A_O2+4.*c.A_CO+5.*c.A_H2O) + ((t.^2-0.298.^2)./2)*(-c.B_But-4.5.*c.B_O2+4.*c.B_CO+5.*c.B_H2O) + ((t.^3-0.298.^3)./3)*(-c.C_But-4.5.*c.C_O2+4.*c.C_CO+5.*c.C_H2O) + ((t.^4-0.298.^4)./4)*(-c.D_But-4.5.*c.D_O2+4.*c.D_CO+5.*c.D_H2O) + ((t.^5-0.298.^5)./5)*(-c.E_But-4.5.*c.E_O2+4.*c.E_CO+5.*c.E_H2O));
dH_Rxn3 = c.H_Rxn3 + ((t-0.298)*(-c.A_But-6.5.*c.A_O2+4.*c.A_CO2+5.*c.A_H2O) + ((t.^2-0.298.^2)./2)*(-c.B_But-6.5.*c.B_O2+4.*c.B_CO2+5.*c.B_H2O) + ((t.^3-0.298.^3)./3)*(-c.C_But-6.5.*c.C_O2+4.*c.C_CO2+5.*c.C_H2O) + ((t.^4-0.298.^4)./4)*(-c.D_But-6.5.*c.D_O2+4.*c.D_CO2+5.*c.D_H2O) + ((t.^5-0.298.^5)./5)*(-c.E_But-6.5.*c.E_O2+4.*c.E_CO2+5.*c.E_H2O));
dH_Rxn4 = c.H_Rxn4 + ((t-0.298)*(-c.A_MA-c.A_O2+4.*c.A_CO+c.A_H2O) + ((t.^2-0.298.^2)./2)*(-c.B_MA-c.B_O2+4.*c.B_CO+c.B_H2O) + ((t.^3-0.298.^3)./3)*(-c.C_MA-c.C_O2+4.*c.C_CO+c.C_H2O) + ((t.^4-0.298.^4)./4)*(-c.D_MA-c.D_O2+4.*c.D_CO+c.D_H2O) + ((t.^5-0.298.^5)./5)*(-c.E_MA-c.E_O2+4.*c.E_CO+c.E_H2O));
dH_Rxn5 = c.H_Rxn5 + ((t-0.298)*(-c.A_MA-3.*c.A_O2+4.*c.A_CO2+c.A_H2O) + ((t.^2-0.298.^2)./2)*(-c.B_MA-3.*c.B_O2+4.*c.B_CO2+c.B_H2O) + ((t.^3-0.298.^3)./3)*(-c.C_MA-3.*c.C_O2+4.*c.C_CO2+c.C_H2O) + ((t.^4-0.298.^4)./4)*(-c.D_MA-3.*c.D_O2+4.*c.D_CO2+c.D_H2O) + ((t.^5-0.298.^5)./5)*(-c.E_MA-3.*c.E_O2+4.*c.E_CO2+c.E_H2O));

%% ODE Equations
dDVdz(1) = K.*(-r.Rxn1 - r.Rxn2 - r.Rxn3);  %delta(n-Butane concentration) over delta(Riser Height). [kmol/s*m]
dDVdz(2) = K.*(r.Rxn1 - r.Rxn4 - r.Rxn5);   %delta(Maleic Anhydride concentration) over delta(Riser Height). [kmol/s*m]
dDVdz(3) = K.*(4.*r.Rxn3 + 4.*r.Rxn5);      %delta(CO2 concentration) over delta(Riser Height). [kmol/s*m]
dDVdz(4) = K.*(4.*r.Rxn2 + 4.*r.Rxn4);      %delta(CO concentration) over delta(Riser Height). [kmol/s*m]
dDVdz(5) = K.*(4.*r.Rxn1 + 5.*r.Rxn2 + 5.*r.Rxn3 + r.Rxn4 + r.Rxn5);    %delta(H2O concentration) over delta(Riser Height). [kmol/s*m] 
dDVdz(6) = K.*(-3.5.*r.Rxn1 - 4.5.*r.Rxn2 - 6.5.*r.Rxn3 - r.Rxn4 - 3.*r.Rxn5);  %delta(O2 concentration) over delta(Riser Height). [kmol/s*m] 
dDVdz(7) = 0;       %delta(inert N2 concentration) over delta(Riser Height). [kmol/s*m]
if s==1     %Adiabatic case. delta(Temperature) over delta(Riser Height). [K/m]
   dDVdz(8) = -K*((dH_Rxn1.*r.Rxn1 + dH_Rxn2.*r.Rxn2 + dH_Rxn3.*r.Rxn3 + dH_Rxn4.*r.Rxn4 + dH_Rxn5.*r.Rxn5)./(nBut.*Cp_But(t) + nMA.*Cp_MA(t) + nCO2.*Cp_CO2(t) + nCO.*Cp_CO(t) + nH2O.*Cp_H2O(t) + nO2.*Cp_O2(t) + nN2.*Cp_N2(t)));
elseif s==2 %Isothermal case (No change in Temperature along riser).
    dDVdz(8)=0;
end
dDVdz(9) = -c.pc.*(1-ep).*c.g;  %delta(Presssure) over delta(Riser Height). [Pa/m]
dDVdz(10) = K*(r.Rxn1);         %delta(Extent of Reaction 1) over delta(Riser Height). [kmol/s*m] 

dDVdz = dDVdz';
end
%%% END FOR ODE SOLVER FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%