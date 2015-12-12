%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single output neural network simulaion
%This is the MAIN function receiving place-cell like input and based on PCA
%like architeture caclulates the output.
%Every network has N1 inputs and NmEC outputs which are independent of
%each-other.
%The main reason for differences between outputs are due to the random
%weights assigned to them initially.
%Yedidyah Dordek, Technion 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
 addpath('functions');
%%  Parameters of simulation
Struct.simdur =5e5; % total simulation time
Struct.Resolution = 75;       %resolution of space, used for determing the velocity
%% Model parameters
% number of grid cells
Struct.NmEC = 1; %for Modules must be more than 1!!!
% % number of place cells 
Struct.N = 25; %number of place cells in one dimension (if not in array this is an avarge)
Struct.N1 = Struct.N^2; % If a non-arrayed input place field is used, its necessary to have suffeicant amount of Place-cells

% learning rate for plasticity from place to grid cells - inital rate. it
% should cool off and lower in time
Struct.epsilon = 1e8;   %the actual learning rate is 1/(t*delta+epsilon)
%cooling parmeter see above.
Struct.delta = 1;
% max weight to a grid cell
Struct.maxWeight = 0.04; %you don't want to give an initial wight with a too high value

%% Sizes of environment
Struct.minx = 0; Struct.maxx =10; %distance in arbitrary units
Struct.miny = 0; Struct.maxy =10;%distance in arbitrary units
%%  Veclocity
Struct.vel  =2.15*Struct.maxx/ Struct.Resolution;
%should reflect the #Place-cells per distance unit. Needs to ~ maxX/(2*N)

% angular velocity. 
% Note: when using inputs of zero mean (like disks or DOGs) use a small angular velocity in order to get a smoother trajectory. 
% When using Gaussian input +direvation, use a large (2*pi) angular velocity for the sake of smoothness in temporal diffrentiation.
Struct.angular = .2;  
%% Type of input
%1. Diff of Gaussians
Struct.placeCellType = 'DOG';

%2. Disks
%Struct.placeCellType = 'Disk';

%3. Gaussians
%Struct.placeCellType = 'Gaussian'; %in this case, need to determine the method of mean zero... (adaptation/diffrentioation)
%% Place-cells properties
% choose the arragnment of place cells. in an array or scattered. NOTE: if
Struct.PcSize = .75; %NOTE!!:when using disks, you need larger sizes. at least 1.5 

% scattered, you'll need more place-cells in a given environment.
Struct.arragmentPC = 'array';
%NOTE: in order for the scattered option to work well (specifically in 
%the case of non restricted weights you need to increase place cell's
%density (and adjust velocity as well...)

% Struct.arragmentPC = 'scattered';

% peak (saturating) value of grid cell activation
Struct.psiSat = 30;
%% MeanZero method (only for option #3 - Gaussians)
Struct.meanZaro = 'adpat'; %adaptation
    % growth rate for onset of adaptation, if used.
    Struct.b1 = 0.5;
    % growth rate of inactivation of adaptation, if used.
    Struct.b2 =Struct.b1/3;

%derivatives
%Struct.meanZaro = 'diff'; %diffrentioations

%% Architecture
%single output
Struct.Arc = 'single';

%multiple output
%Struct.Arc = 'multiple'; %using the symmetric Foldiak's algorithm

%Assymetric architecure
%NOTE: in order for ALL outputs to converge, simulation time needs to be
%long. at least 5e6...
%Struct.Arc = 'sanger'; % using Sanger's algorithm. capable of calc all PCa

%% Output type
%define output fucntion to be sigmoid or linear, NOTE: a large slope (see
%in function) in the sigmoid can cause a quick freeze.

% signoid output
%Struct.output = 'sigmoid';
% Linear output. NOTE: use a rather large slope - help it converge. 
Struct.output = 'linear';
%% Sending all data to function
%activate an interactive plot of the online activity (1 set of temproal weights, avg activity and weights spatial activity)
paint = 1;
%Constrain weights to be nonNegative?
Struct.NonNegativity =1;
% Save input for PCA? 
%NOTE: the longer the simulation -> more data needed.
% with 16GB  of memory duration of simulation should not exceed 4e5, with 100 outputs, 900 inputs.
Struct.saveInput = 1;
if Struct.saveInput
%Allocate memory to speed up process
  Struct.TemporalInput = zeros(Struct.N1,Struct.simdur);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%send the data for simulation !!!!
[Struct] = PCANetworkFunc(Struct, paint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PCA and Spatial Maps
[Struct] = PCAEvecCalc(  Struct,Struct.TemporalInput,Struct.NmEC);
%clear  Struct.TemporalInput  %clear some memory if needed...


%calculate the Maps, Gridness, angles, etc.
%weightsMap holds the var requested to calc.

weightsMap=Struct.J; %J are weights of the network
%weightsMap =Struct.PCANNEvec; %calc for NonNegative PCA
%weightsMap =Struct.PCAEvec; %calc for "regular" PCA

%in order to plot a map you need to send in 
 [Struct.outputMap, Struct.autocorrMap, Struct.Gridness60,...
     Struct.Gridness90, Struct.minAngles, Struct.GridSpacing] = Maps(Struct,weightsMap);

 %% Plot results & Autocorrelations
numberOfPlots=2;
n_mat=mulMatCross(zeros(size(Struct.outputMap,1)),zeros (size(Struct.outputMap,1)));
plotResults( Struct.outputMap,n_mat,numberOfPlots)
 
 
 
 