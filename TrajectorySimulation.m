clear all;
close all; 

%% Simulation parameters

% Number of Ag atoms to simulate
N = 1e7;

% Magnetic field status (on=1/off=0)
Bon = 1; 

% Plot during simulation
ploton = 0;

% Save final workspace and plots
saveresults = 1;

%% Physical constants and experimental parameters

% Bohr magneton [J/s]
mu_b = 9.27*10^(-24); 

% Mass of Ag [kg]
m = 1.7912* 10^(-25);

% Oven apperature radius [m]
oven_inlet_radius =  5.6419e-04;
% y-position of the oven [m] (from paper (1.75+3.3+2.5)/100 )
oven_y = -.0765; 
% S1 apperature radius [m]
rslit1_inlet =  3.09e-05; 
% y-position of S1 [m] 
slit1y = -5.15e-2;
% S2 apperature width [m]
slt_w = 8e-4;
% S2 apperature height [m]
slt_h = 3.5e-5;
% y-position of S2 [m]
slit2y = -1.8e-2;
% Center of O, S1, and S2 cutplanes: z_c in paper [m]
slt_z =-1.9e-4 - 30e-6;

%% Kinematic details of atoms starting at oven

% Initialize dynamic variables for coordinates, velocities and accelarations
r = zeros(N,3);     % Coordinates [m]
v = zeros(N,3);     % Velocities [m/s]
a = zeros(N,3);     % Acceletations [m/s^2]

% Set initial coordinates to oven slit [m]
r(:,2) = oven_y;
r(:,1) = oven_inlet_radius/3.*randn(N,1);
r(:,3) = oven_inlet_radius/3.*randn(N,1) + slt_z;
rinit = r;

% Angular Distribution deflecting towards the first slit
% The atoms that would not pass throught the first slit are excluded
% Initialize the direction of the velocity of each atom
alpha = rslit1_inlet/abs(2.5e-2); %Cone vertex angle
theta = (2*rand(N,1)-1)* alpha/2;
phi = rand(N,1)* 2*pi; 
u = [0 slit1y slt_z] - r;
u = u./vecnorm(u')';
u(:,[1 2 3]) = u(:,[1 3 2]); 
ind = u(:,3) > 1-eps;
u_prime(ind,1) = sin(theta(ind)).*cos(phi(ind));
u_prime(ind,2) = sin(theta(ind)).*sin(phi(ind));
u_prime(ind,3) = sign(u(ind,3)).*cos(theta(ind));
ind = ~ind;
u_prime(ind,1) = sin(theta(ind)).*(u(ind,1).*u(ind,3).*cos(phi(ind)) - u(ind,2).*sin(phi(ind)) )./sqrt(1-u(ind,3).^2) + u(ind,1).*cos(theta(ind));
u_prime(ind,2) = sin(theta(ind)).*(u(ind,2).*u(ind,3).*cos(phi(ind)) + u(ind,1).*sin(phi(ind)) )./sqrt(1-u(ind,3).^2) + u(ind,2).*cos(theta(ind));
u_prime(ind,3) = -sqrt(1-u(ind,3).^2).*sin(theta(ind)).*cos(phi(ind)) + u(ind,3).*cos(theta(ind));
u_prime(:,[1 2 3]) = u_prime(:,[1 3 2]);

% Initialize speed distribution
vstart = 625;                               % [m/s] minimum speed
vend = 750;                                 % [m/s] minimum speed
vabs = (vend-vstart)* rand(N,1) + vstart;   % Uniform sample
vinit = vabs .* u_prime;                    % Initial velocities
v = vinit;                                  % Dynamic variable

% Spin (random, +1 or -1)
spin = randi(2,N,1)*2-3; 
% spin = rand(N,1)*2-1; 
% Alive status (atom's that have hit an object will be set to 0)
alive = ones(N,1);
% Arrival status (atom's that have hit the detector will be set to 1)
finished = zeros(N,1);

%%
tic

% Load Magnetic field data from an external file. 
% Simulated in Comsol Multiphysics and exported. 
% Data file is too large for online storage. Please reach out for details.
comsoldatagridded = importfile('BfieldGrad.txt');

z_SG = 30e-6*(1/sin(32*pi/180)-1);

% Copy the grid values entered in comsol to here
fieldmeshx = [-0.5:0.01:0.5]'/100;          % [m]  
fieldmeshy = [-5.1:0.01:5]'/100;            % [m]  
fieldmeshz = [-0.5:0.01:0.5]'/100;          % [m]
meshsize = [length(fieldmeshx) length(fieldmeshy) length(fieldmeshz)]; 

% Reshape the field data into 3D grid
x_comsol = reshape(comsoldatagridded(:,1),meshsize)/100;     % [T] 
y_comsol = reshape(comsoldatagridded(:,2),meshsize)/100;     % [T]
z_comsol = reshape(comsoldatagridded(:,3),meshsize)/100;     % [T] 
Bx = reshape(comsoldatagridded(:,4),meshsize);     % [T] 
By = reshape(comsoldatagridded(:,5),meshsize);     % [T]
Bz = reshape(comsoldatagridded(:,6),meshsize);     % [T] 
% Reshape the field gradient data into 3D grid
dBzdx = reshape(comsoldatagridded(:,13),meshsize); % [T/m]
dBzdy = reshape(comsoldatagridded(:,14),meshsize); % [T/m]
dBzdz = reshape(comsoldatagridded(:,15),meshsize); % [T/m]
% The following terms are irrelevant due to z axis quantization
%dBxdx = reshape(comsoldatagridded(:,7),meshsize);
%dBxdy = reshape(comsoldatagridded(:,8),meshsize);
dBxdz = reshape(comsoldatagridded(:,9),meshsize);
%dBydx = reshape(comsoldatagridded(:,10),meshsize);
%dBydy = reshape(comsoldatagridded(:,11),meshsize);
dBydz = reshape(comsoldatagridded(:,12),meshsize);

Bx_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,Bx,'linear','none');
By_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,By,'linear','none');
Bz_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,Bz,'linear','none');
dBzdx_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,dBzdx,'linear','none');
dBzdy_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,dBzdy,'linear','none');
dBzdz_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,dBzdz,'linear','none');
dBxdz_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,dBxdz,'linear','none');
dBydz_interpolant = griddedInterpolant(x_comsol,y_comsol,z_comsol,dBydz,'linear','none');


%% Main simulation

toc 

% Final time to stop simulation [s]
tf = 1e-3;
% Temporal step size [s]
dt = 1.5133e-6/2;
ts = 0:dt:tf;

if saveresults
    strnow = string(datetime('now','TimeZone','local','Format','y-M-d_HH-mm-ss_z'));
    mkdir(strnow)
    copyfile([mfilename('fullpath'),'.m'],strnow+"/executedscript.m");
end

reverseStr = '';

% Set temporary dynamic variables
rtemp = r;
for it = 1:length(ts)
    t = ts(it);
    
    % Set all accelerations to zero
    a(:,1:3) = 0;
    % If magnetic field ON, calculate accelaration for alive atoms
    if Bon == 1
        
        dBxdz_atom = dBxdz_interpolant(rtemp(alive==1,1),rtemp(alive==1,2),rtemp(alive==1,3));
        dBydz_atom = dBydz_interpolant(rtemp(alive==1,1),rtemp(alive==1,2),rtemp(alive==1,3));
        dBzdz_atom = dBzdz_interpolant(rtemp(alive==1,1),rtemp(alive==1,2),rtemp(alive==1,3));
        dBxdz_atom(isnan(dBxdz_atom)) = 0;
        dBydz_atom(isnan(dBydz_atom)) = 0;
        dBzdz_atom(isnan(dBzdz_atom)) = 0;
        
        % Equation of motion
        a(alive==1,3) = mu_b*spin(alive==1).*dBzdz_atom/m;
        a(alive==1,2) = mu_b*spin(alive==1).*dBydz_atom/m;
        a(alive==1,1) = mu_b*spin(alive==1).*dBxdz_atom/m;

        % If outside the magnet region, set a = 0 
        a(~(rtemp(:,2)> -.0175 & rtemp(:,2)<.0175),:) = 0;

    end
    
    % Update velocity
    v = v + a*dt;
    
    % Update next projected coordinates
    rtemp = r + v*dt;

    % 'Kill' the atom if it 'hits' an object
    % First slit conditional (kill atoms that interact)
    alive(rtemp(:,2)>slit1y & r(:,2)<slit1y & ...
       (rtemp(:,1)).^2+(rtemp(:,3)-slt_z).^2>rslit1_inlet^2 ...
        ) = 0;
    % Second slit conditional
    alive(rtemp(:,2)>slit2y & r(:,2)<slit2y & ...
        (p1(rtemp(:,1),slt_h,slt_w,slt_z)<rtemp(:,3 ) | ...
        p2(rtemp(:,1),slt_h,slt_w,slt_z)>rtemp(:,3)) ...
        ) = 0;

    % 'Finish' the atoms that 'hit' the Detector plate
    finished(rtemp(:,2)> .0175) = 1;
    
    % Update the next coordinate for atoms that are 'alive' and 'unfinished'
    r(alive==1 & ~finished,:) = rtemp(alive==1 & ~finished,:);
    
    % Termination condition when less than 5% are left to finish
    if size(finished(~finished& alive==1))/size(finished(finished& alive==1))<0.05
        alive(~finished& alive==1) = 0; % Kill unfinished
        fprintf(['\n']);
        disp("Less than 5% left running");
        break;
    end
    
    msg = sprintf('t = %f us', t*1e6);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    if ploton
        % Live Propagation projected in 3D
        if ~exist("fig",'var')
            fig = figure(101);
            figs1 = scatter3(r(alive==1,1),r(alive==1,2),r(alive==1,3),5,'k');
            hold on;
            figs2 = scatter3(r(alive==0,1),r(alive==0,2),r(alive==0,3),1,'r');
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            hold off;
            title("t=  "+ string(t));
            xlim([-6e-4, 6e-4]);
            ylim([-8e-2, 2e-2]);
            drawnow;
        else
            figs1.XData = r(alive==1,1);
            figs1.YData = r(alive==1,2);
            figs1.ZData = r(alive==1,3);
            figs2.XData = r(alive==0,1);
            figs2.YData = r(alive==0,2);
            figs2.ZData = r(alive==0,3);
            fig.CurrentAxes.Title.String = "t =  "+ num2str(t,3) + ' s';
            drawnow;
        end
    end
end

%% Plot and save results

figure;
scatter(r(alive==1,1).*1000,r(alive==1,3).*1000,.5,'k');
xlabel('x (mm)');
ylabel('z (mm)');
grid on;
axis equal;
box on;
drawnow

if saveresults
    saveas(gcf,strnow+"\fig1.png")
end

figure;
scatter(r(alive==1,1).*1000,r(alive==1,3).*1000,.5,v(alive==1,2));
colormap jet;
xlabel('x (mm)');
ylabel('z (mm)');
grid on;
axis equal;
ylim([-0.45 0.1])
xlim([-0.9 0.9])
box on;
clbr = colorbar;
clbr.Label.String = 'v_y (m/s)';
clbr.Label.FontSize = 12;
drawnow

if saveresults
    saveas(gcf,strnow+"\fig2.png")
    
    aalive = a(alive==1,:);
    ralive = r(alive==1,:);
    valive = v(alive==1,:);
    spinalive = spin(alive==1,:);

    save(strnow+"\workspace" + ".mat",'aalive','ralive','valive','spinalive')
end


toc

%% Helper functions

% Parametric functions describing Slit S2
function p1 = p1(x,a,b,c)
    p1 = -a/2/(b/2)^2 .* x.^2 + c +a/2;
end

function p2 = p2(x,a,b,c)
    p2 = a/2/(b/2)^2 .* x.^2 + c -a/2;
end 


%Importing function
function ComsolFieldExport = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  COMSOLFIELDEXPORT = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the numeric data.
%
%  COMSOLFIELDEXPORT = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  ComsolFieldExport = importfile("ComsolFieldExport.txt", [10, Inf]);
%

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [10, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 28);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["VarName1", "Version", "COMSOL", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"];
opts.SelectedVariableNames = ["VarName1", "Version", "COMSOL", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["VarName1", "Version", "COMSOL", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15"], "ThousandsSeparator", ",");

% Import the data
ComsolFieldExport = readtable(filename, opts);

%% Convert to output type
ComsolFieldExport = table2array(ComsolFieldExport);
end


