%% loads data for fitting

%% Data ([ATP] = 8, 4, 2 mM)
ATP_c = [8, 4, 2]; % mM for each column
Data_ATP = [0		56.4048	63.6074	61.3192
        0.5		51.812	51.8626	47.4794
        1		37.4459	35.9182	31.387
        2		17.8025	13.5516	10.2112
        3		11.443	8.34	6.3895
        4		6.2643	2.8669	2.5781
        5		3.2759	1.6526	1.594
        6		2.212	1.1823	1.2117];

% Anthony Baker 2021 06 15 isovelocity fit Filip.xlsx
Data_ATP_AB2021 = [...
0	67.39416105	65.42069119
0.5	60.98847686	43.94730968
1	42.70587933	27.84044045
2	18.85391508	8.062961483
3	12.59560277	5.10199899
4	6.442581818	2.930592431
5	3.513565502	1.870419908
6	1.85557286	1.041556855];


% Non-zero velocities
vel = (-Data_ATP(:,1)); % ML/s

%% plot the datasets together
% plot(Data_ATP(:, [2 3 4]), Data_ATP(:, 1), 'o-', Data_ATP_AB2021(:, [2 3]), Data_ATP_AB2021(:, 1), 'o--', LineWidth=2, MarkerSize=8);hold off;
% legend('8mM','4mM', '2mM','8mM AB2021','2mM AB2021')
% plot(Data_ATP(:, [2 4]), Data_ATP(:, 1), 'o-', Data_ATP_AB2021(:, [2 3]), Data_ATP_AB2021(:, 1), 'o--', LineWidth=2, MarkerSize=8);hold off;
% legend('8mM','2mM','8mM AB2021','2mM AB2021')

%% Fmax (normalized.) versus [MgATP] (mM) from Ebus et al.(2001)
% iso_data = ...
%     [0.0098736     0.019874      0.04959     0.098478      0.49024        5.063
%       1.5925       1.6826       1.5898       1.4657       1.2884      0.99732];

iso_data = ...
    [0.01     0.02      0.05     0.1      0.5        5.0
      1.5925       1.6826       1.5898       1.4657       1.2884      0.99732];
  
MgATP_iso = iso_data(1,:);

F_data = iso_data(2,:).*57;
F_data_r = iso_data(2,:); % relative data

%% Ktr data from Beard et al. - 8, 4, 2 mM
MgATP = [8 4 2];
Ktr_mean = [37.7928 29.0 25.8033];
Ktr_err  = [1.9308  1.30    2.0167];

%% Km range
% Yamashita et al 1994: ADP Inhibits the Sliding Velocit of Fluorescent Actin Filaments on Cardiac and Skeletal Myosins
% they got a value of 0.043
% HOwever the in vitro preparation underestimates the Km. Thus, we
% extrapolate the error to form a range

k_m_ADP0 = [0.04 0.1];

% extrapolation for nonzero MgADP
% k_m_=0.043(1 + [MgADP]/194)

