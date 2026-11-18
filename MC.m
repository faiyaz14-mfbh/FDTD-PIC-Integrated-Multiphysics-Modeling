
%close all;

  clear;
  clc;
  
 %%
 
 % TIME TO REACH THE END
 15E-6/1E5
  
%% 
  eps0 = 8.854e-12;
  qe = 1.602e-19;
  k = 1.381e-23;
  amu = 1.661e-27;
  m = 32.0 * amu;
  
  %% excitation freq at the beg of the 2DEG
  
f=1e11;

%%
v_drift     = 1.15e5; %m/s
len_x_2DEG  = 20e-6;

N_2DEG          = 2000;
dx_2DEG       = len_x_2DEG/N_2DEG;
X             = (1:N_2DEG)*dx_2DEG;


periodic_file_dump_on=1;
filename   = ['OUtput_'];


%% new particle

 
% initializations 
    part_max  =  1.7e5;
    np_insert =  10;
    np          =  0;
    tt_insert =  1000;  % TT_max is about 2 M
    dt           =  1.1793e-17;

    
     np_i      = np_insert/tt_insert;  % inserted number of particles per iteration

%     % scale time for faster MC with same MF
      dt_MC     = 1000*dt;
      np_i_MC  = 1000*np_i;  % make sure this is integer
    
    
    
          
    X_2DEG = dx_2DEG*(1:N_2DEG+1);
    vel_times_PICPC           = zeros(1, N_2DEG);
    J_e                              = zeros(1, N_2DEG);
    
    % 2DEG density
    
    n_sh_cm2 = 1e12;
    n_sh_m2  = 6e8*1e4;
    

    t_2DEG  =  10e-9;
    
    % Macro factor calculation
    PIC                    =np_i_MC*dx_2DEG/(v_drift*dt_MC);
    macro_sheet_density    = PIC / (dx_2DEG);  % m^-2
    MF                     = n_sh_m2 /macro_sheet_density;  
    
    
%%
    
  part_x = zeros(part_max,1);
  part_v = zeros(part_max,1);
 
  e               = zeros(1, N_2DEG);
  PIC_per_cell    = zeros(1, N_2DEG+1);



  fprintf ( 1, '\n' );
  fprintf ( 1, '    Time  Particles\n' );
  fprintf ( 1, '\n' );

 step_num  = 1600*11;
for tt = 1: step_num

% =========  1. CALCULATE THE CHARGE DENSITY.  
PIC_per_cell    = zeros(1, N_2DEG+1);
vel_times_PICPC    = zeros(1, N_2DEG+1);

    for p = 1 : np

      fi = 1.0 + part_x(p) / dx_2DEG;
      i  = floor ( fi );
      hx = fi - i;

      PIC_per_cell(i)     = PIC_per_cell(i)       + (1.0 - hx );
      vel_times_PICPC(i)     = vel_times_PICPC(i)       + (1.0-hx)*part_v(p);
      
      PIC_per_cell(i+1)   = PIC_per_cell(i+1)     +      hx;
      vel_times_PICPC(i+1)   = vel_times_PICPC(i+1)     +      hx*part_v(p);

    end

%====== 4. GENERATE NEW PARTICLES

        if (part_max - np <= np_insert)
             np_insert = part_max - np;  % random velocity and n matrices may have issue with this
        end

        % particle position initialization in the very first cell.
            part_x(np+1:np+np_insert,1) = rand(np_insert, 1 ) * dx_2DEG;    
        % normal velocty distribution initialization
            part_v(np+1:np+np_insert,1) = normrnd(v_drift, v_drift/200, [np_insert, 1]); 
            np = np + np_insert;

% point excitation sinosoidal.
%     e      = zeros(1, N_HD);
%     f=1e11;
%     e(5)   = 1e9*sin(2*pi*f*dt*step);
   
% random velocity at half grind points using normal dstribution

  p = 1;
% update each particle's position and velocity
    while ( p <= np )
      
      xx       = 1.0 + part_x(p) / dx_2DEG;
      part_pos = floor ( xx );
      
% E field and Particle position arrays
%      E   E   E   E   E   E   E   E
%    n   n   n   n   n   n   n   n   n
% therefoe partcle at ith location interacts with e(i)
      
      f = -qe * e(part_pos);%* cos(2*pi*6e10*dt*step_num);
      acc = f / m;
      part_v(p) = part_v(p) + acc * dt_MC;
      part_x(p) = part_x(p) + part_v(p) * dt_MC;
        
% Absorption:
% Kill the particle by replacing it with the last particle.
      if ( part_x(p) >= len_x_2DEG)
        part_x(p)= part_x(np);
        part_v(p) = part_v(np);
       disp('deleted');
        np = np - 1;
        p = p - 1;
      end
      p = p + 1;
      
    end
np;


n_sh_macro  = vel_times_PICPC/dx_2DEG ;
J2_e       =  -n_sh_macro(1:N_2DEG)*MF*qe/t_2DEG;

if(mod(tt, 1000)==0)
        subplot(3,1,1); plot(X, PIC_per_cell(1:N_2DEG), '--b'); title('No. of particles'); drawnow;
        subplot(3,1,2); plot(part_x, part_v, 'b');    title('position versus velocity'); drawnow;
        subplot(3,1,3); plot(X, -J2_e, '--m');         title('J'); drawnow;
end
        if(periodic_file_dump_on==1)
           if mod(tt,1000)==0
               filNameWorkspace = [filename num2str(tt) '.mat'];
               save(filNameWorkspace);
           end
        end
    
end


%%
figure;
subplot(2,1,1); plot(X, PIC_per_cell(1:N_2DEG), 'om');
subplot(2,1,2); scatter(part_x, part_v, 'o');
figure;
subplot(2,1,1); plot(vel_times_PICPC)


%%
% 
% 
% [f,xi] = ksdensity(part_v);
% figure; plot(xi,f); 
% 
% 
% [f,xi] = ksdensity(part_x);
% figure; plot(xi,f); 
% 

