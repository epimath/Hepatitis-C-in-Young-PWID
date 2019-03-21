%% PWID Syringe Sharing Sharing Contact Matrices (Smith MK et al. STEP Study, Baltimore)

% Full citation of the Smith MK et al. STEP study:
% Smith MK, Graham M, Latkin CA, Mehta SH, Cummings DAT. (2018).
% Quantifying potentionally infectious sharing patterns among people who
% inject drugs in Baltimore, USA. Epidemiology and Infection, 146,
% 1845-1853.

% Contact matrix used to fit model in:
% Gicquelais RE, Foxman B, Coyle J, Eisenberg MC.(2019).
% Hepatitis C transmission in young people who inject drugs: insights using
% a dynamic model informed by state public health surveillance. Epidemics.
% https://doi.org/10.1016/j.epidem.2019.02.003

% This script re-creates the matrix of deviations from the proportional
% mixing assumption based on a study of PWID in Baltimore, MD. Each number
% in this matrix represents the magnitude of syringe sharing observed in a
% study of PWID in Baltimore relative to expected sharing that would be 
% expected under proportional mixing (i.e. reflects assortative
% syringe sharing by age).

% These matrices are used to model syringe sharing in a hepatitis C 
% transmission model in Gicquelais et al. 

% Age groups (and matrix rows/columns) include:
% <25	25-29	30-34	35-39	40-44   45-49   50-54   55+;

% A second matrix of 1s represents the proportional mixing assumption.

%% Magnitude of assortative over proportional mixing
assort = [2.68  6.31    3.02    0.63    0.7     0       0.72    0;
        0.49    2.17    1.3     0.35    0.14    0.34    1.29    0.28;
        0.23    1.19    1.22    1.83    1.14    0.78    1.54    0.92;
        0.06    0.53    0.9     2.16    1.24    1.35    1.56    1.8;
        0.09    0.3     0.48    1.55    1.88    2.68    2.72    0.79;
        0       0.16    0.35    0.77    0.99    2.51    2.5     1.54;
        0.01    0.07    0.1     0.46    0.38    1.29    2.52    2.55;
        0.03    0.09    0.17    0.7     0.26    1.67    3.24    7.74];

% Combine age groups as necessary for age structure in Gicquelais et al.

%calc median for <25 and 25-29 sharing with 30+
young_shr_old=[assort(3:end,1),assort(3:end,1),assort(3:end,2)];
young_shr_old1=median(young_shr_old,1);

old_shr_young=[assort(1,3:end);assort(1,3:end);assort(2,3:end)];
old_shr_young1=median(old_shr_young,2);

old_shr_old=assort(3:end,3:end);
old_shr_old1=median(old_shr_old(:));

avgc=[assort(1,1) assort(1,1) assort(1,2) old_shr_young1(1,1);
            assort(1,1) assort(1,1) assort(1,2) old_shr_young1(2,1);
            assort(2,1) assort(2,1) assort(2,2) old_shr_young1(3,1);
            young_shr_old1(1,1) young_shr_old1(1,2) young_shr_old1(1,3) old_shr_old1(1,1)];

%% Create proportional mixing matrix (ones)
prop=ones(4,4);

%% Create minimum deviance from prop mixing matrix for LHS sampling
minc=[min(assort(1,1),prop(1,1)) min(assort(1,1),prop(1,2)) min(assort(1,2),prop(1,3)) min(min(old_shr_young(1,:)),prop(1,4));
      min(assort(1,1),prop(2,1)) min(assort(1,1),prop(2,2)) min(assort(1,2),prop(2,3)) min(min(old_shr_young(2,:)),prop(2,4));
      min(assort(2,1),prop(3,1)) min(assort(2,1),prop(3,2)) min(assort(2,2),prop(3,3)) min(min(old_shr_young(3,:)),prop(3,4));
      min(min(young_shr_old(:,1)),prop(4,1)) min(min(young_shr_old(:,2)),prop(4,2)) min(min(young_shr_old(:,3)),prop(4,3)) min(min(old_shr_old(:)),prop(4,4))];

%% Create maximum deviance from prop mixing matrix for LHS sampling
maxc=[max(assort(1,1),prop(1,1)) max(assort(1,1),prop(1,2)) max(assort(1,2),prop(1,3)) max(max(old_shr_young(1,:)),prop(1,4));
      max(assort(1,1),prop(2,1)) max(assort(1,1),prop(2,2)) max(assort(1,2),prop(2,3)) max(max(old_shr_young(2,:)),prop(2,4));
      max(assort(2,1),prop(3,1)) max(assort(2,1),prop(3,2)) max(assort(2,2),prop(3,3)) max(max(old_shr_young(3,:)),prop(3,4));
      max(max(young_shr_old(:,1)),prop(4,1)) max(max(young_shr_old(:,2)),prop(4,2)) max(max(young_shr_old(:,3)),prop(4,3)) max(max(old_shr_old(:)),prop(4,4))];


c=avgc;

%%
clearvars('avgc','assort','old_shr_old','old_shr_old1','old_shr_young','old_shr_young1','prop','young_shr_old','young_shr_old1')

