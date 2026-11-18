%==========================================================================
% FDTD_WG_PIC.m
%
% 2D FDTD solver (TEz polarization: Ex, Ey, Hz) coupled to a 1D
% particle-in-cell (PIC / Monte Carlo) model of the electrons in a 2DEG
% channel.
%
% STRUCTURE
%   A parallel-plate waveguide (PEC walls at top and bottom), with PML on
%   the left and right ends. Layer stack, top to bottom:
%     PEC wall | air | gate layer with metal grating | barrier | 2DEG |
%     substrate (fine mesh) | substrate (coarse mesh) | PEC wall
%
% COUPLING (every time step)
%   1. FDTD updates Ex, Ey, Hz (Yee scheme, split field Hz = Hzx + Hzy).
%   2. Ex along the 2DEG row is interpolated to each particle and pushes
%      it with F = -q*E.
%   3. The particles' AC velocities are deposited onto the grid
%      (cloud-in-cell) to give the channel current density J.
%   4. J is fed back into the Ex update (Ampere's law), so the 2DEG
%      current radiates into the grid.
%
% REQUIRES
%   MC_Current_1p0e5_dx_10nm_L_5um_delV_200.mat  (initial particle state
%   and 2DEG parameters; see the "Particle initialization" section)
%
% OUTPUTS
%   *.bin                       probe fields and particle data (raw doubles)
%   Output_<timestamp><step>.mat periodic workspace snapshots
%   Output_<timestamp>.mat       final workspace
%
% ORIGINAL NOTES
%   This file is for testing only the current in a 2DEG using an EM-MC
%   coupled model. The test simply confirms whether the required current
%   density is observed in the channel. There is no grating, source or
%   drain used for this problem.
%   NOTE: the current version DOES include a 29-strip grating and an Ex
%   source in the grating gaps (see "Geometry" and "Source and
%   measurement probes"), so this note seems out of date.
%
%   2DEG thickness increased to 10 nm for a faster simulation.
%
%   MF (macro factor) definition revised based on the injection model,
%   i.e. on conservation of particles upon injection.
%==========================================================================

clear;
% Reference frequency [Hz], used below for lambda and omega.
% NOTE: the name f is reused later. The source frequency is reset to
% 300 GHz inside the time loop, and f is also used for the particle force
% (-qe*E) in the particle mover.
f =3e12;
close all;

%% Output file naming
% Build a unique suffix from the current date/time (e.g. _20261117213605)
% so each run writes new output files instead of overwriting old ones.
timestamp = fix(clock);
uniq_identifier = '_';
for i=1:length(timestamp)
    uniq_identifier = [uniq_identifier num2str(timestamp(i))];
end

filename = ['Output_' uniq_identifier];
%% Simulation mode switches
% recording_on          : 1 = write probe fields / particle data to .bin files
% mesh_only             : 1 = build and plot the mesh/geometry, skip the time loop
% periodic_file_dump_on : 1 = save the whole workspace every period_file_dump
% period_file_dump      : snapshot interval [s]

recording_on   = 1;
mesh_only      = 0;

periodic_file_dump_on    = 1;
period_file_dump         =  0.5e-12;


%% FDTD init
% Physical constants (SI units)

epsi_o    =  8.854e-12;       % F/m
mu_o      =  4*pi*1e-7;
c         =  1/sqrt(mu_o*epsi_o);
lambda    =  (c/f);                 % m N_HD  1 THz
omega     =  2*pi*c/lambda;
eta       =  sqrt(mu_o/epsi_o);

% Relative permittivities of the substrate and barrier (9.5, roughly GaN)
epsi_sub  = 9.5;
epsi_barr = 9.5;


%%  non-uniform mesh

% ---- y-directed (vertical) mesh, listed top to bottom ----
% For each layer: len_y_* = thickness [m], dy_* = cell size [m],
%                 M_*     = number of cells in that layer.
%   metal  : top PEC wall
%   air2   : air gap above the gate
%   gate   : gate / grating layer
%   barr   : barrier (AlGaN)
%   2DEG   : a single cell row holding the 2DEG (dy = its thickness)
%   subs1  : substrate, fine mesh near the channel
%   subs2  : substrate, coarser mesh
%   metal2 : bottom PEC wall

len_y_metal=1e-6;  dy_metal=100e-9;           M_metal=round(len_y_metal/dy_metal);
len_y_air2=1e-6;   dy_air2=50e-9;                 M_air2=round(len_y_air2/dy_air2);
len_y_gate=30e-9;  dy_gate=10e-9;               M_gate=round(len_y_gate/dy_gate);
len_y_barr=30e-9;  dy_barr=10e-9;                M_barr=round(len_y_barr/dy_barr);
len_y_2DEG=10e-9;  dy_2DEG=len_y_2DEG;    M_2DEG=round(len_y_2DEG/dy_2DEG);        
len_y_subs1=1e-6;  dy_subs1=10e-9;            M_subs1=round(len_y_subs1/dy_subs1);              
len_y_subs2=1e-6;  dy_subs2= 50e-9;           M_subs2=round(len_y_subs2/dy_subs2);
len_y_metal2=1e-6; dy_metal2=100e-9;        M_metal2=round(len_y_metal2/dy_metal2);


% 2DEG thickness [m]; converts sheet current [A/m] to current density [A/m^2]
t_2DEG   = dy_2DEG;

% Total number of cells in y
M        =  M_metal+M_air2+M_gate+M_barr+M_2DEG+M_subs1+M_subs2+M_metal2;

% ---- x-directed (horizontal) mesh, listed left to right ----
%   PML | coarse | transition | 2DEG region (fine, 10 nm) | transition |
%   coarse | PML
% For each region: len_x_* = length [m], dx_* = cell size [m],
%                  N_*     = number of cells.
% len_PML is used only for the (disabled) y-PML; len_PML_x is the x-PML.
len_PML     = 80e-6;       
len_PML_x   = 80e-6;

len_x_1=len_PML_x;    dx_1=10e-6;      N_1=round(len_x_1/dx_1);
len_x_2=80e-6;         dx_2=5e-6;       N_2=round(len_x_2/dx_2);
len_x_3=1e-6;         dx_3=100e-9;     N_3=round(len_x_3/dx_3);
len_x_2DEG=5e-6;      dx_2DEG=10e-9;   N_2DEG=round(len_x_2DEG/dx_2DEG);
len_x_4=1e-6;         dx_4=100e-9;     N_4=round(len_x_4/dx_4);
len_x_5=80e-6;         dx_5=5e-6;       N_5=round(len_x_5/dx_5);
len_x_6=len_PML_x;    dx_6=10e-6;      N_6=round(len_x_6/dx_6);

% Total number of cells in x
N            =          N_1+N_2+N_3+N_2DEG+N_4+N_5+N_6;

%% Cell-size maps
% DX(i,j), DY(i,j) = width / height of the cell at row i (y), column j (x).
% Stored as full M-by-N matrices so the update equations can use them
% element by element.
DX      = cat(2,ones(M,N_1)*dx_1,        ...
                ones(M,N_2)*dx_2,        ...
                ones(M,N_3)*dx_3,        ...
                ones(M,N_2DEG)*dx_2DEG,  ...
                ones(M,N_4)*dx_4,        ... 
                ones(M,N_5)*dx_5,        ...
                ones(M,N_6)*dx_6); 

DY      = cat(1,   ...
                ones(M_metal,N)*dy_metal,   ...
                ones(M_air2,N)*dy_air2,   ...
                ones(M_gate,N)*dy_gate,   ...
                ones(M_barr,N)*dy_barr,   ...
                ones(M_2DEG,N)*dy_2DEG,   ...
                ones(M_subs1,N)*dy_subs1, ...
                ones(M_subs2,N)*dy_subs2, ...
                ones(M_metal2,N)*dy_metal2   ...
                );
            
% Time step from the 2D Courant (CFL) stability limit, using the smallest
% cells in x and y.
dt      = 1/(c*sqrt(1/min(DX(1,1:N))^2 + 1/min(DY(1:M,1))^2));

% Physical coordinates of the cells (for plotting)
    X           =  cumsum(DX(1,1:N));
    Y           =  cumsum(DY(1:M,1));
    [XX YY]     =  meshgrid(X,Y);

    % diag_figure: a map used only to plot the geometry. Each region gets a
    % tag value: 1 = background, 8 = PML, 30 = substrate, 50 = barrier,
    % 60 = 2DEG, 100 = PEC walls, 200 = grating / probes,
    % 190 = input/output probes, +/-180 = source.
    diag_figure =  ones(M,N);

%% material properties
% Start with free space everywhere; dielectrics are added in "Geometry".
mu       =  ones(M,N)*mu_o;
epsi     =  ones(M,N)*epsi_o;


%% Particle (MC / PIC) initialization
% Load the initial particle state and 2DEG parameters. Variables used
% below that come from this .mat file include:
%   part_x, part_v      particle positions [m] and velocities [m/s]
%   np, part_max        number of active particles, size of the arrays
%   np_i, np_insert,    injection settings
%   tt_insert
%   v_drift             drift velocity [m/s]
%   n_sh_m2             physical sheet electron density [m^-2]
load('MC_Current_1p0e5_dx_10nm_L_5um_delV_200.mat');
% Electron charge [C], free electron mass [kg], effective mass (0.2*m0, GaN)
   qe = 1.602e-19;
   m_o = 9.1093e-31;
   me  = 0.2*m_o;
 
% x-positions of the 2DEG cell edges, from the start of the channel [m]
    X_2DEG_edge                   = dx_2DEG*(0:N_2DEG);
    
% 2DEG cell centers [m]; the very last one is an extra point past the end
     X_2DEG_CellCenter              = dx_2DEG*(0:N_2DEG)+dx_2DEG/2;
     
% Arrays on the N_2DEG+1 edge points:
%   j_sh_density_macro      macro-particle current at the nodes
%   j_sh_density_macro_1by2 same, interpolated to the cell edges
%   e                       Ex sampled along the 2DEG [V/m]
%   J                       physical current density fed to FDTD [A/m^2]
     j_sh_density_macro_1by2  = zeros(1, N_2DEG+1);
     j_sh_density_macro       = zeros(1, N_2DEG+1);
     e                        = zeros(1, N_2DEG+1);
     J                        = zeros(1, N_2DEG+1);   

    % Split each particle's velocity into a DC part (its injection velocity)
    % and an AC part (the change caused by the field). Only the AC part is
    % deposited as current, so the steady drift current does not radiate.
    % part_E holds the field seen by each particle.
    part_v_dc  =  part_v;
    part_v_ac  =  zeros(size(part_v));
    part_E     =  zeros(part_max,1);
    
    % Macro-particle scaling:
    %   PIC = np_i*dx_2DEG/(v_drift*dt): expected macro-particles per cell,
    %         i.e. injection rate times the transit time across one cell.
    %   macro_sheet_density = macro-particles per unit length of channel.
    %   MF  = macro factor: real electrons represented by one macro-particle,
    %         chosen so the macro-particle density matches n_sh_m2.
    PIC                    = np_i*dx_2DEG/(v_drift*dt);
    macro_sheet_density    = PIC / (dx_2DEG);  % m^-2
    MF                     = n_sh_m2 /macro_sheet_density;

%% Initialize field arrays (all fields start at zero)
Ex    = zeros(M,N);
Ey    = zeros(M,N);
Hzy   = zeros(M,N);
% Hz is split into Hzy and Hzx (Berenger split-field PML): Hz = Hzy + Hzx
Hzx   = zeros(M,N);
Hz    = zeros(M,N);


% Electric conductivities: used for the PML, and to model PEC metal as a
% very high conductivity.
sigma_ey  = zeros(M,N);
sigma_ex  = zeros(M,N);

% PML thicknesses in cells. Only PML_x is used; PML_y and PML_y_t belong
% to the commented-out y-PML below.
PML_x   =  round(len_PML_x/dx_1);
PML_y   =  round(len_PML/dy_metal);
PML_y_t =  round(len_y_subs2/dy_subs2 + len_y_barr/dy_barr + (len_PML-len_y_subs2-len_y_barr)/dy_metal);

%% PML conditions
% Electrical conductivity.
% Anisotropic conductivity (sigma_ex, sigma_ey) is used to build the PML
% at the borders. The same arrays are also used in the main region, where
% sigma_ex = sigma_ey and sigma_mx = sigma_my.
%
% x-PML: conductivity grows as (depth/PML_x)^m, up to 100 S/m, toward the
% left and right ends. The top and bottom are closed by PEC walls instead
% (see "Geometry"), so the y-PML is commented out.
 m=2;
% Grading profile: polynomial of order m, maximum value 100
gradPMLx= ((1:PML_x)' / PML_x).^m * 100;
% Fill the right-hand PML (columns near N) and mirror it on the left.
% sigma_ex(:,N-PML_x) is still zero here, so (... + 1)*gradPMLx is just
% gradPMLx.
for u=N-PML_x+1:N
    sigma_ex(:,u-1)        = (sigma_ex(:,N-PML_x) +1) * (gradPMLx(u-(N-PML_x)));
    sigma_ex(:,N-u+1)      =  (sigma_ex(:,N-PML_x) +1) * (gradPMLx(u-(N-PML_x)));
    diag_figure(:,u-1)     = 8*ones(size(diag_figure(:,u)));
    diag_figure(:,N-u+1)   = 8*ones(size(diag_figure(:,N-u+1)));
end

% y-PML (disabled: top and bottom are PEC walls)
% m=4;
% gradPMLy= ((1:PML_y)' / PML_y).^m * 10000;
% for u=M-PML_y+1:M
% %   sigma_ey(u,X_out_side_HEMT)      = (sigma_ey(M-PML_y,X_out_side_HEMT) +1) * (gradPMLy(u-(M-PML_y)));
%      sigma_ey(u,:)      = (sigma_ey(M-PML_y,:) +1) * (gradPMLy(u-(M-PML_y)));
%      sigma_ey(M-u+1,:)  = (sigma_ey(M-PML_y,:) +1) * (gradPMLy(u-(M-PML_y)));
%      diag_figure(u,:)       = diag_figure(u,:)  + 8*ones(size(diag_figure(u,:)));
%      diag_figure(M-u+1,:)   = diag_figure(M-u+1,:) +8*ones(size(diag_figure(M-u+1,:)));
% end

% gradPMLy_t= ((1:PML_y_t)' / PML_y_t).^m * 1000;
% for u=M-PML_y_t+1:M
%     sigma_ey(u,X_out_side_HEMT)   =  (sigma_ey(M-PML_y_t,X_out_side_HEMT) +1) * (gradPMLy_t(u-(M-PML_y_t)));
% end

% Magnetic conductivities (same indexing as Hz)
sigma_mx   = zeros(M,N); % analogous to Hz
sigma_my   = zeros(M,N); % analogous to Hz

% PML matching condition: sigma_m = sigma_e*mu/epsi makes the PML
% impedance match the interior, so waves enter it without reflecting.
%sigma_my(M-PML_y:M,:)  =  sigma_ey(M-PML_y:M,:) .* mu(M-PML_y:M,:)./epsi(M-PML_y:M,:);
%sigma_my(1:PML_y,:)    =  sigma_ey(1:PML_y,:)   .* mu(1:PML_y,:)./epsi(1:PML_y,:);
sigma_mx(:,N-PML_x:N)  =  sigma_ex(:,N-PML_x:N) .* mu(:,N-PML_x:N)./epsi(:,N-PML_x:N);
sigma_mx(:,1:PML_x)    =  sigma_ex(:,1:PML_x)   .* mu(:,1:PML_x)./epsi(:,1:PML_x);

% Coefficients for a 1D wave used as the SF/TF reference (not used now).
% These should be calculated on the bare grid with PML only, without any
% features (conductivity) introduced.


%% Geometry

% 2DEG: one cell row directly below the barrier.
% POS1_2deg_x = [first last] column of the channel, POS1_2deg_y = its row.
    POS1_2deg_x   = [N_1+N_2+N_3   N_1+N_2+N_3+N_2DEG];
    POS1_2deg_y   = M_metal+M_air2+M_gate+M_barr+1;
    diag_figure(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2)) = 60*ones(size(diag_figure(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2))));
    
% Waveguide walls: PEC modeled as a very high conductivity (1e20 S/m) in
% the top and bottom metal layers.
    sigma_ey(1:M_metal,:)         = 1*1e20;                   % pec TOP
    sigma_ex(1:M_metal,:)         = 1*1e20;                   % pec TOP 
    sigma_ey(M-M_metal2+1:M, :)   = 1*1e20;                   % pec BOTTOM 
    sigma_ex(M-M_metal2+1:M, :)   = 1*1e20;                   % pec BOTTOM 
    
    diag_figure(1:M_metal,:)      = 100*ones(size(diag_figure(1:M_metal,:)));
    diag_figure(M-M_metal2+1:M,:) = 100*ones(size(diag_figure(M-M_metal2+1:M, :)));
    
% Grating: num_of_grat PEC strips, each gate_w wide, separated by
% gate_gap, in the gate layer. The first strip starts 6 cells into the
% channel. Each loop pass draws one strip, then steps to the next one.
    gate_y_start  = M_metal+M_air2+1;
    gate_y_end    = M_metal+M_air2+M_gate;

    grating_x_start= POS1_2deg_x(1)+6;
    num_of_grat  =  29;
    gate_gap     =  80e-9;
    gate_w       =  90e-9;
    grating_x_start_tmp = grating_x_start;

    for gr=1:num_of_grat
        sigma_ey(gate_y_start:gate_y_end, grating_x_start_tmp:grating_x_start_tmp+round(gate_w/dx_2DEG))   = 1*1e20;        % pec strip 
        sigma_ex(gate_y_start:gate_y_end, grating_x_start_tmp:grating_x_start_tmp+round(gate_w/dx_2DEG))   = 1*1e20;  
        diag_figure(gate_y_start:gate_y_end, grating_x_start_tmp:grating_x_start_tmp+round(gate_w/dx_2DEG)) = 200*ones(size(diag_figure(gate_y_start:gate_y_end,grating_x_start_tmp:grating_x_start_tmp+round(gate_w/dx_2DEG))));   % end
        grating_x_start_tmp = grating_x_start_tmp + (gate_gap+gate_w)/dx_2DEG;
    end
    


% Dielectric substrate (2DEG row + both substrate layers), only under the
% channel region in x
    sub_y_start  = M_metal+M_air2+M_gate+M_barr+1;  sub_y_end   = sub_y_start+M_2DEG+M_subs1+M_subs2-1;
    epsi(sub_y_start:sub_y_end,POS1_2deg_x(1):POS1_2deg_x(2))        = epsi(sub_y_start:sub_y_end,POS1_2deg_x(1):POS1_2deg_x(2))*epsi_sub;
    diag_figure(sub_y_start:sub_y_end,POS1_2deg_x(1):POS1_2deg_x(2)) = 30*ones(size(diag_figure(sub_y_start:sub_y_end,POS1_2deg_x(1):POS1_2deg_x(2))));

% Dielectric barrier, only under the channel region in x.
% NOTE: this uses epsi_sub, not epsi_barr. Both are 9.5 now so there is no
% difference, but there will be if epsi_barr is changed.
    barr_y_start = M_metal+M_air2+M_gate+1;         barr_y_end  = M_metal+M_air2+M_gate+M_barr;
    epsi(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2))        = epsi(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2))*epsi_sub;
    diag_figure(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2)) = 50*ones(size(diag_figure(barr_y_start:barr_y_end,POS1_2deg_x(1):POS1_2deg_x(2))));

%% Source and measurement probes
    
% Source: Ex is driven in each gap between grating strips, with the sign
% alternating from gap to gap (+1, -1, +1, ...). source_profile stores
% this +/-1 pattern, indexed from grating_x_start.


%     source_profile = 1-2*(sigma_ey(gate_y_start+1,POS1_2deg_x(1):POS1_2deg_x(2))...
%         /max(sigma_ey(gate_y_start+1,POS1_2deg_x(1):POS1_2deg_x(2))))

    grating_x_start_tmp = grating_x_start;
    for gr=1:num_of_grat-1    
        if(mod(gr,2)==0)
          source_profile(grating_x_start_tmp-grating_x_start+round(gate_w/dx_2DEG) ...
                        :grating_x_start_tmp-grating_x_start+round(gate_w/dx_2DEG)+round(gate_gap/dx_2DEG))=-1;
                    
        else
         source_profile(grating_x_start_tmp-grating_x_start+round(gate_w/dx_2DEG) ...
                        :grating_x_start_tmp-grating_x_start+round(gate_w/dx_2DEG)+round(gate_gap/dx_2DEG))=1;
                      
        end
        % (No semicolon on the next line, so MATLAB prints it every iteration.)
        grating_x_start_tmp = grating_x_start_tmp + (gate_gap+gate_w)/dx_2DEG
    end
    

    % Source location: the columns covered by source_profile, at mid-height
    % of the gate layer
    source_x = grating_x_start+1:grating_x_start+size(source_profile,2); % N_1+3;
    source_y =  M_metal+M_air2+round(M_gate/2);
    diag_figure(source_y, source_x) = 180*source_profile;
   
    % Probe lines where fields are recorded:
    %   meas_IP   vertical, just left of the channel (gate layer to substrate)
    %   meas_OP   vertical, just right of the channel (same height range)
    %   meas_OP2  vertical, near the right PML, full waveguide height
    %   meas_barr horizontal, through the middle of the barrier, along the channel
    meas_IP_x = POS1_2deg_x(1)-3;
    meas_IP_y = gate_y_start:sub_y_end;
    diag_figure(meas_IP_y, meas_IP_x) = 190;
    
    meas_OP_x = POS1_2deg_x(2)+3;
    meas_OP_y = gate_y_start:sub_y_end;
    diag_figure(meas_OP_y, meas_OP_x) = 190;
 
    meas_OP2_x = N-N_6-4;
    meas_OP2_y = M_metal+2:M-M_metal2-2;
    diag_figure(meas_OP2_y, meas_OP2_x) = 200;

    meas_barr_x = POS1_2deg_x(1):POS1_2deg_x(2);
    meas_barr_y = M_metal+M_air2+M_gate+round(M_barr/2);
    diag_figure(meas_barr_y, meas_barr_x) = 200;

%% Coefficients for the update equations
% Standard lossy-FDTD coefficients (sigma*E averaged over the time step).
% For each field: c_*e* multiplies the old field, c_*h* the curl term.
% Ex (updated from dHz/dy) uses sigma_ey; Ey (from dHz/dx) uses sigma_ex.
% E-field update coefficients
c_exe    = (1-sigma_ey*dt./(2*epsi)) ./ (1+sigma_ey*dt./(2*epsi));
c_exh    = 1./(1+sigma_ey*dt./(2*epsi));
c_eye    = (1-sigma_ex*dt./(2*epsi)) ./ (1+sigma_ex*dt./(2*epsi));
c_eyh    = 1./(1+sigma_ex*dt./(2*epsi));

% H-field update coefficients (separate sets for the Hzy and Hzx parts)
% c_hzh    = (1-sigma_m*dt./(2*mu)) ./ (1+sigma_m*dt./(2*mu));
c_hzhy    = (1-sigma_my*dt./(2*mu)) ./ (1+sigma_my*dt./(2*mu));
c_hzey    = 1./(1+sigma_my*dt./(2*mu));
c_hzhx    = (1-sigma_mx*dt./(2*mu)) ./ (1+sigma_mx*dt./(2*mu));
c_hzex    = 1./(1+sigma_mx*dt./(2*mu));

%% grid: illustration

%  Ey  Hz  Ey  Hz  Ey  Hz  Ey  Hz
%      Ex      Ex      Ex      Ex
%  Ey  Hz  Ey  Hz  Ey  Hz  Ey  Hz
%      Ex      Ex      Ex      Ex
%  Ey  Hz  Ey  Hz  Ey  Hz  Ey  Hz
%      Ex      Ex      Ex      Ex

% ^  this is y direction
% >  this is x direction

% E and H are updated 0.5*dt apart (leapfrog); successive H updates are
% dt apart.
% In every field matrix, rows correspond to y and columns to x.

%% In this code we are dealing with the TEz polarized wave:
%  Hz, Ex, Ey components are being considered.




%% Geometry plot (cell-index axes)
figure;
h= pcolor(diag_figure);
set(h, 'EdgeColor', 'none');
% set(gca, 'XLim',[N_1 N-N_1], 'YLim',[M_1 M-M_1]);
colorbar;

%% Geometry plot (physical x-y axes, m)

        figure;
        h= pcolor(XX, YY,diag_figure);
        set(h, 'EdgeColor', 'none');
        
%      set(gca, 'XLim',[len_x_1 X(end)-len_x_5], 'YLim',[len_y_1 Y(end)-len_y_6]);
        colorbar;

%% Old test settings (kept for reference)

%  sigma_ex(M_1:M_1+M_2-1, N_1:N_1+N_2-1) = 1e20;
%  sigma_ey(POS_2deg_y,POS_2deg_x(1):POS_2deg_x(2)-1) = 1e10/2;
%  sigma_ex(PEC_y_start:M-2,PML_x+10:N-PML_x-11)      = 1e10/5;
%  sigma_ex(POS_2deg_y,POS_2deg_x(1):POS_2deg_x(2)-1) = 1e10/2;
%  figure;    pcolor(XX,YY,sigma_ex);
   
    
% T_toal and k are not used anywhere below.
T_toal = 20000;
k      = 10;

%% Open binary files for probe / particle dumps
% Each file name gets this run's unique suffix. Data are appended as raw
% doubles every tt_rec steps (see "Runtime plots" in the time loop).

if (recording_on)

        file_Ex_IP = ['Ex_IP_' uniq_identifier '.bin'];
        file_Ey_IP = ['Ey_IP_' uniq_identifier '.bin'];
        file_Hz_IP = ['Hz_IP_' uniq_identifier '.bin'];
        
        file_Ex_OP = ['Ex_OP_' uniq_identifier '.bin'];
        file_Ey_OP = ['Ey_OP_' uniq_identifier '.bin'];
        file_Hz_OP = ['Hz_OP_' uniq_identifier '.bin'];
            
        file_Ex_OP2 = ['Ex_OP2_' uniq_identifier '.bin'];
        file_Ey_OP2 = ['Ey_OP2_' uniq_identifier '.bin'];
        file_Hz_OP2 = ['Hz_OP2_' uniq_identifier '.bin'];
        
        file_Exo_barr = ['Ex_barr' uniq_identifier '.bin'];
        file_Eyo_barr = ['Ey_barr' uniq_identifier '.bin'];
        file_Hzo_barr = ['Hz_barr' uniq_identifier '.bin'];

        file_e        = ['e_' uniq_identifier '.bin'];
        file_n        = ['n_' uniq_identifier '.bin'];
        file_J        = ['J_' uniq_identifier '.bin'];
        file_part_x   = ['Part_x_' uniq_identifier '.bin'];
        file_part_v   = ['Part_v_' uniq_identifier '.bin'];
      
  % open the files
        fExi1  = fopen(file_Ex_IP,'w');
        fEyi1  = fopen(file_Ey_IP,'w');
        fHzi1  = fopen(file_Hz_IP,'w');
                    
        fExo1  = fopen(file_Ex_OP,'w');
        fEyo1  = fopen(file_Ey_OP,'w');
        fHzo1  = fopen(file_Hz_OP,'w');
        
        fExo2  = fopen(file_Ex_OP2,'w');
        fEyo2  = fopen(file_Ey_OP2,'w');
        fHzo2  = fopen(file_Hz_OP2,'w');
        
        fExo_barr  = fopen(file_Exo_barr,'w');
        fEyo_barr  = fopen(file_Eyo_barr,'w');
        fHzo_barr  = fopen(file_Hzo_barr,'w');
        
        fe     = fopen(file_e,'w');    
        fn     = fopen(file_n,'w'); 
        fJ     = fopen(file_J,'w'); 
        fPart_x     = fopen(file_part_x,'w');    
        fPart_v     = fopen(file_part_v,'w');         
end


%% Time iteration loop
% t_total : total simulated time [s] (100e-10 s = 10 ns)
% tt_max  : number of FDTD time steps
% tt_rec  : plot and record every tt_rec steps
% (t_record_step is computed but not used.)
t_total = 100e-10;
tt_max  = fix(t_total/dt);
tic

tt_rec  = 2000;
t_record_step = dt*tt_rec;

if (~mesh_only)
figure;
for tt =1:tt_max  % tt = 1:150 K
    
%% SF/TF boundary conditions for plane wave source.
% 1D plane wave for SF/TF  (hypothetical vertical propagating wave)
%      x1                  = ((tt*dt-tao)/sigma1).*((tt*dt-tao)/sigma1);
%      Hzy_1D(1)           = exp(0.5.*(-x1))/377;   % guassian pulse with max E-field as 1 V/m
%      Hzy_1D(1)           = A*pi*fc*(tt*dt-tao).*exp(-2*(pi*fc*(tt*dt-tao)).^2)/377;   % guassian pulse with max E-field as 1 V/m     

     % Only this first sample of the 1D SF/TF source is still computed; it is
     % not used, since the SF/TF corrections are commented out.
     Hzy_1D(1)           = sin(2*pi*f*(tt)*dt)/377;

%     Ex_1D(1:sr_len-1)   = c_exe_1D(1:sr_len-1)'.*Ex_1D(1:sr_len-1) ...
%                         + c_exh_1D(1:sr_len-1)' .* (dt./dy_1D(1:sr_len-1)'./epsi_o).*(Hzy_1D(2:sr_len)-Hzy_1D(1:sr_len-1));
%                     
%     Hzy_1D(2:sr_len)    = c_hzhy_1D(2:sr_len)'.*Hzy_1D(2:sr_len) ...
%                             + c_hzey_1D(2:sr_len)'.*(dt./dy_1D(2:sr_len)'./mu_o).*(Ex_1D(2:sr_len)-Ex_1D(1:sr_len-1));

      
% ---------------- 2D FDTD update ----------------
% Source: set Ex in the grating gaps to a 300 GHz sine (amplitude 100 V/m)
% with the alternating +/-1 profile, before the regular Ex update.

   f=300e9;
   Ex(source_y, source_x)    =  1e2*sin(2*pi*f*dt*tt)*source_profile ;
 
   % Ex update from dHz/dy (Ampere's law, y-propagating part)
   Ex(1:M-1,1:N)  =    c_exe(1:M-1,1:N).*Ex(1:M-1,1:N)  ...
                      + c_exh(1:M-1,1:N).*(dt./DY(1:M-1,1:N)./epsi(1:M-1,1:N)).*(Hz(2:M,1:N) - Hz(1:M-1,1:N));
                  
% Ey update from dHz/dx (Ampere's law, x-propagating part)

   Ey(1:M,2:N)    =    c_eye(1:M,2:N).*Ey(1:M,2:N)   ...
                     - c_eyh(1:M,2:N).*(dt./DX(1:M,2:N)./epsi(1:M,2:N)).*(Hz(1:M,2:N) - Hz(1:M,1:N-1));

% SF/TF correction                 
%    Ex(y_sr,1:N)   =   Ex(y_sr,1:N)  + ...
%                                  c_exh(y_sr,1:N).*   ...
%                                  (dt./DY(y_sr,1:N)./epsi(y_sr,1:N)).*Hzy_1D(1);                       

% Coupling to the 2DEG: subtract the channel current J from Ex along the
% 2DEG row (Ampere's law: dE/dt = curl(H)/eps - J/eps). This is how the
% particle current radiates into the FDTD grid.

Ex(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2)) = Ex(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2)) ...
                                             - J*dt./epsi(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2));
                                 
                                                         
% Hz update                                      
% Hzy update from dEx/dy (Faraday's law, y-propagating part)
  Hzy(2:M,1:N-1)  = c_hzhy(2:M,1:N-1).*Hzy(2:M,1:N-1)    ...
                    + c_hzey(2:M,1:N-1).*(dt./DY(2:M,1:N-1)./mu(2:M,1:N-1)).*(Ex(2:M,1:N-1) - Ex(1:M-1,1:N-1));
                
% Hzx update from dEy/dx (Faraday's law, x-propagating part)
  Hzx(2:M,1:N-1)   =  c_hzhx(2:M,1:N-1).*Hzx(2:M,1:N-1)    ...
                     - c_hzex(2:M,1:N-1).*(dt./DX(2:M,1:N-1)./mu(2:M,1:N-1)).*(Ey(2:M,2:N) - Ey(2:M,1:N-1));
                
  % Copy column 1 into the last column of Hzy/Hzx. (Hz(:,N) itself is
  % overwritten by the boundary condition below.)
  Hzy(2:M,N)   =  Hzy(2:M,1);
  Hzx(2:M,N)   =  Hzx(2:M,1);
  
% SF/TF correction for Hzy
%    Hzy(y_sr,x_sr(1:end-1))  =  Hzy(y_sr,x_sr(1:end-1))  ... 
%                                +   c_hzey(y_sr,x_sr(1:end-1)).*(dt./DY(y_sr,x_sr(1:end-1))./mu(y_sr,x_sr(1:end-1))).*Ex_1D(1);

%   Hzy(y_sr,1:N-1)     =  Hzy(y_sr,1:N-1)  ...
%                                +   c_hzey(y_sr,1:N-1).*(dt./DY(y_sr,1:N-1)./mu(y_sr,1:N-1)).*Ex_1D(1);
                           
% Total Hz = Hzy + Hzx
   Hz(2:M,1:N-1)    =   Hzy(2:M,1:N-1) +  Hzx(2:M,1:N-1);
   
% Neumann boundary on Hz (dHz/dn = 0, equivalent to PEC) on the top row
% and the last column
   Hz(1,1:N-1)   = Hz(2,1:N-1);
   Hz(1:M,N)     = Hz(1:M,N-1);



%% Field probes (old version, disabled)
%    Ex_input(tt, :)        =   Ex(y_sr+1,x_meas1);
%    Hz_input(tt, :)        =   Hz(y_sr+1,x_meas1);
%    Ex_output(tt, :)       =   Ex(y_meas1,x_meas1);
%    Hz_output(tt, :)       =   Hz(y_meas1,x_meas1);
%    Ex_output2(tt, :)      =   Ex(y_meas2,x_meas2);
%    Hz_output2(tt, :)      =   Hz(y_meas2,x_meas2);
   

%% Runtime plots and data recording (every tt_rec steps)
% Panels: Ey map, macro-particles per cell, Ex along the 2DEG, current.
    if mod(tt,tt_rec)==0
          subplot(4,1,1); pcolor(Ey);  caxis([-100 100]*10); drawnow; colorbar;
          subplot(4,1,2); plot(X_2DEG_edge, PIC_per_cell, '-R'); title('No. of particles'); drawnow;         
          subplot(4,1,3); plot(e);%, vel_times_PICPC(1:N_2DEG)./PIC_per_cell(1:N_2DEG), '-B');    title('vel'); drawnow;
          subplot(4,1,4); plot(j_sh_density_macro, '-B');         title('J'); drawnow;

      fprintf('FDTD sim Time = %f ps;  Sim. time %f s\n', tt*dt/1e-12, toc);

% Append the current probe values and particle data to the .bin files
        if (recording_on)

            fwrite(fExi1,Ex(meas_IP_y,meas_IP_x), 'double');
            fwrite(fEyi1,Ey(meas_IP_y,meas_IP_x), 'double');
            fwrite(fHzi1,Hz(meas_IP_y,meas_IP_x), 'double');

            fwrite(fExo1,Ex(meas_OP_y,meas_OP_x), 'double');
            fwrite(fEyo1,Ey(meas_OP_y,meas_OP_x), 'double');
            fwrite(fHzo1,Hz(meas_OP_y,meas_OP_x), 'double');  
        
            fwrite(fExo2,Ex(meas_OP2_y,meas_OP2_x), 'double');
            fwrite(fEyo2,Ey(meas_OP2_y,meas_OP2_x), 'double');
            fwrite(fHzo2,Hz(meas_OP2_y,meas_OP2_x), 'double');    

            fwrite(fExo_barr, Ex(meas_barr_y,meas_barr_x), 'double');
            fwrite(fEyo_barr, Ey(meas_barr_y,meas_barr_x), 'double');
            fwrite(fHzo_barr, Hz(meas_barr_y,meas_barr_x), 'double');
        
             fwrite(fe,e, 'double');    
             fwrite(fn,PIC_per_cell, 'double'); 
             fwrite(fJ,J, 'double'); 
             fwrite(fPart_x, part_x, 'double');    
             fwrite(fPart_v, part_v, 'double');
                    
        end
     
   end
    

%% 2DEG channel: particle (PIC) update

% time iteration for the MC mover
  
%====== 1.  field coupled to channel
% Sample Ex at the N_2DEG+1 cell edges along the 2DEG row

e(1:N_2DEG+1) = Ex(POS1_2deg_y,POS1_2deg_x(1):POS1_2deg_x(2));

%====== 2. GENERATE NEW PARTICLES at every tt_insert

  if mod(tt,tt_insert)==0
        % Inject np_insert new particles every tt_insert steps, capped so the
        % total never exceeds part_max.
        if (part_max - np <= np_insert)
             np_insert = part_max - np;  % random velocity and n matrices may have issue with this
        end
% Place the new particles at random positions in the first cell
            part_x(np+1:np+np_insert,1) = rand(np_insert, 1 ) * dx_2DEG;    
% Velocities: Gaussian around v_drift with std v_drift/200; this is also
% stored as their DC velocity
            part_v(np+1:np+np_insert,1) = normrnd(v_drift, v_drift/200, [np_insert, 1]); 
            part_v_dc(np+1:np+np_insert,1) = part_v(np+1:np+np_insert,1);
            np = np + np_insert;
            disp('inserted')
  end

  
% Grid layout along the channel:
%   e : N_2DEG+1 values, Ex at the cell edges (first edge to last edge)
%   n : N_2DEG+1 values, particles per cell; the last cell is extra, so
%       the first N_2DEG values are the real cells


% Linearly interpolate Ex from the edges to each particle's position
part_E = interp1(X_2DEG_edge, e(1:N_2DEG+1), part_x);
%e_center = interp1(X_2DEG_edge, e(1:N_2DEG+1), X_2DEG_CellCenter);

%====== 3. update each particle's position and velocity
p = 1;
    
% A while loop is used because particles that leave the channel are
% removed in place, so np changes during the loop.
while ( p <= np )
       % Cell index of the particle (only used by the commented-out e_center line)
       xx       = 1.0 + part_x(p) / dx_2DEG;
       part_pos = floor ( xx );

% Push: force -> acceleration -> velocity -> position.
% part_v_ac accumulates only the field-driven velocity change.
% NOTE: this reuses the variable name f (the frequency) for the force.

      f    = -qe * part_E(p);
      %f         = -qe *e_center(part_pos);
      acc       = f / me;
      part_v(p)    = part_v(p) + acc * dt;
      part_v_ac(p) = part_v_ac(p)+acc * dt;
      part_x(p)    = part_x(p) + part_v(p) * dt;
        
% Absorption:
% Kill the particle by replacing it with the last particle.
     if (part_x(p) >= len_x_2DEG)||(part_x(p) <= 0)
        part_x(p)= part_x(np);
        part_v(p) = part_v(np);
        part_v_dc(p) = part_v_dc(np);
        part_v_ac(p) = part_v_ac(np);
        part_E(p)    =part_E(np);
        np = np - 1;
        p = p - 1;
     end
      p = p + 1;
      
end
    
% =========  4. CALCULATE THE CURRENT DENSITY.  

% Deposit particles onto the grid with cloud-in-cell (linear) weighting:
% a particle at fractional position hx inside cell i gives weight (1-hx)
% to node i and hx to node i+1. Both the particle count (PIC_per_cell)
% and the AC-velocity sum (vel_times_PICPC) are deposited.
PIC_per_cell                = zeros(1, N_2DEG+1);
vel_times_PICPC             = zeros(1, N_2DEG+1);

    for p = 1 : np
        
      fi = 1.0 + part_x(p) / dx_2DEG;
      i  = floor ( fi );
      hx = fi - i;

      % particles in each cell
      PIC_per_cell(i)     = PIC_per_cell(i)   + (1.0 - hx );
      vel_times_PICPC(i)  = vel_times_PICPC(i)  ...     
                               + (1.0-hx)*(part_v_ac(p));
      
      PIC_per_cell(i+1)      = PIC_per_cell(i+1)  +  hx;
      vel_times_PICPC(i+1)   = vel_times_PICPC(i+1)     ...
                              +  hx*(part_v_ac(p));
     
    end

% There are N_2DEG+1 cells, we ignore the first and the last cell, 
% so we have N_2DEG-1 cells
    
% Macro-particle current per unit length at the nodes, then interpolated
% from cell centers to the interior edges; the two end edges copy their
% neighbors.
 j_sh_density_macro                = vel_times_PICPC/dx_2DEG;
 j_sh_density_macro_1by2(2:N_2DEG) =  ...
                                    interp1(X_2DEG_CellCenter(1:N_2DEG),...
                                            j_sh_density_macro(1:N_2DEG), ...
                                            X_2DEG_edge(2:N_2DEG));
 
 j_sh_density_macro_1by2(1)        = j_sh_density_macro_1by2(2);
 j_sh_density_macro_1by2(N_2DEG+1) = j_sh_density_macro_1by2(N_2DEG);
  
% Old averaging method (kept for reference):
% % Ignore 1st and last cell, now we have N_2DEG-1 cells with N_2DEG edges
% % current at theses N_2DEG edges is
% % j_sh_density_macro_1by2(2:N_2DEG-1) ...
% %    = 0.5 * (j_sh_density_macro(2:end-2) + j_sh_density_macro(3:end-1));
% % 

%%j_sh_density_macro_1by2(N_2DEG) = j_sh_density_macro_1by2(N_2DEG-1);
 

 % Convert to physical current density [A/m^2]:
 %   * MF     -> real electrons per macro-particle
 %   * qe     -> charge (the minus sign is for the negative electron charge)
 %   / t_2DEG -> spread the sheet current over the 2DEG thickness
 J        =  -j_sh_density_macro_1by2*MF*qe/t_2DEG;


%% Periodic workspace snapshot (every period_file_dump seconds)

        if(periodic_file_dump_on==1)
           if mod(tt-1,round(period_file_dump/dt))==0
               filNameWorkspace = [filename num2str(tt) '.mat'];
               save(filNameWorkspace);
           end
        end
end
end

%% Save the final workspace and close the dump files
        save(filename);
        if (recording_on)             
            
             fclose(fExi1);
             fclose(fEyi1);
             fclose(fHzi1);
                    
             fclose(fExo1);
             fclose(fEyo1);
             fclose(fHzo1);
             
             fclose(fExo2);
             fclose(fEyo2);
             fclose(fHzo2);

             fclose(fExo_barr);
             fclose(fEyo_barr);
             fclose(fHzo_barr);

             fclose(fe);    
             fclose(fn); 
             fclose(fJ); 
             fclose(fPart_x);    
             fclose(fPart_v);
        end


disp('end of calculations')


    