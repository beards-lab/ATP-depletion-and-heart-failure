%% loads data for fitting

%% Data ([ATP] = 8, 4, 6 mM)
Data_ATP = [0		56.4048	63.6074	61.3192
        0.5		51.812	51.8626	47.4794
        1		37.4459	35.9182	31.387
        2		17.8025	13.5516	10.2112
        3		11.443	8.34	6.3895
        4		6.2643	2.8669	2.5781
        5		3.2759	1.6526	1.594
        6		2.212	1.1823	1.2117];

Data_ATP = [0		56.4048	63.6074	61.3192
        1		37.4459	35.9182	31.387
        2		17.8025	13.5516	10.2112
        6		2.212	1.1823	1.2117];

ML = 1.1; % half sarcomere length (microns)
% Non-zero velocities
vel = (-Data_ATP(:,1)).*ML; % micron per sec

%% Fmax (normalized.) versus [MgATP] (mM) from Ebus et al.(2001)
iso_data = ...
    [0.0098736     0.019874      0.04959     0.098478      0.49024        5.063
      1.5925       1.6826       1.5898       1.4657       1.2884      0.99732];
MgATP_iso = iso_data(1,:);

F_data = iso_data(2,:).*57;

%% Ktr data from Beard et al. - 8, 4, 2 mM
MgATP = [8 4 2];
Ktr_mean = [37.7928 29.0 25.8033];
Ktr_err  = [1.9308  1.30    2.0167];

%% Set up environment
MgADP = 0; 
Pi    = 0; 
