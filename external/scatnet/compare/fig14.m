% This code reproduces the results of the scattering transform for fig12 in
% the paper "PHASE SPACE SKETCHING FOR CRYSTAL IMAGE ANALYSIS BASED ON 
% SYNCHROSQUEEZED TRANSFORMS" by Jianfeng Lu and H, aizhao Yang. 
%
% By Haizhao Yang, 2018


%% load the database
clear; close all;
% NOTE : the following line must be modified with the path to the
% uiuc database in YOUR system.
path_to_db = './scatnet/data/crystal/fig14/';
src = syn_src(path_to_db);
db_name = 'syn14';

use_precomputed_scattering = 0; % change to 0 to skip computation of scattering

grid_train = [36]; % number of training for classification
nb_split = 1; % number of split for classification




%% ---------------------------------------------------
%% ----------------- trans_scatt ---------------------
%% ---------------------------------------------------






%% compute scattering of all images in the database
tic;
feature_name = 'trans_scatt';
precomputed_path = sprintf('./scatnet/precomputed/%s/%s.mat', db_name, feature_name);
if (use_precomputed_scattering)
    load(precomputed_path);
else
    %configure scattering
    options.J = 5; % number of octaves
    options.Q = 1; % number of scales per octave
    options.M = 2; % scattering orders
    
    % build the wavelet transform operators for scattering
    Wop = wavelet_factory_2d_pyramid();
    
    % a function handle that
    %   - read the image
    %   - resize it to 200x200
    %   - compute its scattering
    fun = @(filename)(scat(imresize_notoolbox(imreadBW(filename)...
        ,[200 200]), Wop));
    
    % compute all scattering
    % (800 seconds on a 2.4 Ghz Intel Core i7)
    trans_scatt_all = srcfun(fun, src);
    
    % a function handle that
    %   - format the scattering in a 3d matrix
    %   - remove margins
    %   - average accross position
    % (10 seconds on a 2.4 Ghz Intel Core i7)
    fun = @(Sx)(mean(mean(remove_margin(format_scat(Sx),0),2),3));
    trans_scatt = cellfun_monitor(fun ,trans_scatt_all);
    
    % save scattering
    save(precomputed_path, 'trans_scatt'); 
end
% format the database of feature
db = cellsrc2db(trans_scatt, src);
timeScat1 = toc;

%% classification
tic;
rsds_classif2(db, db_name, feature_name, grid_train, nb_split);
timeCls1 = toc;



%% ---------------------------------------------------
%% --------------- roto_trans_scatt ------------------
%% ---------------------------------------------------





%% compute scattering of all images in the database
tic;
feature_name = 'roto_trans_scatt';
precomputed_path = sprintf('./scatnet/precomputed/%s/%s.mat', db_name, feature_name);
if (use_precomputed_scattering)
    load(precomputed_path);
else
    % configure scattering
    options.J = 5; % number of octaves
    options.Q = 1; % number of scales per octave
    options.M = 2; % scattering orders
    
    % build the wavelet transform operators for scattering
    Wop = wavelet_factory_3d_pyramid();
    
    % a function handle that
    %   - read the image
    %   - resize it to 200x200
    %   - compute its scattering
    fun = @(filename)(scat(imresize_notoolbox(imreadBW(filename),[200 200]), Wop));
    % (1000 seconds on a 2.4 Ghz Intel Core i7)
    roto_trans_scatt_all = srcfun(fun, src);
    
    % a function handle that
    %   - format the scattering in a 3d matrix
    %   - remove margins
    %   - average accross position
    fun = @(Sx)(mean(mean(remove_margin(format_scat(Sx),1),2),3));
    % (10 seconds on a 2.4 Ghz Intel Core i7)
    roto_trans_scatt = cellfun_monitor(fun ,roto_trans_scatt_all);
    
    % save scattering
    save(precomputed_path, 'roto_trans_scatt');
end
% format the database of feature
db = cellsrc2db(roto_trans_scatt, src);
timeScat2 = toc;

%% classification
tic;
rsds_classif2(db, db_name, feature_name, grid_train, nb_split);
timeCls2 = toc;



%% ---------------------------------------------------
%% ------------- roto_trans_scatt_log ----------------
%% ---------------------------------------------------



tic;
feature_name = 'roto_trans_scatt_log';
precomputed_path = sprintf('./scatnet/precomputed/%s/%s.mat', db_name, feature_name);
if (use_precomputed_scattering)
    load(precomputed_path);
else
    % a function handle that
    %   - format the scattering in a 3d matrix
    %   - take the logarithm
    %   - remove margins
    %   - average accross position
    fun = @(Sx)(mean(mean(log(remove_margin(format_scat(Sx),1)),2),3));
    roto_trans_scatt_log = cellfun_monitor(fun ,roto_trans_scatt_all);
    
    %save scattering
    save(precomputed_path);
end
% format the database of feature
db = cellsrc2db(roto_trans_scatt_log, src);
timeScat3 = toc;

%% classification
tic;
rsds_classif2(db, db_name, feature_name, grid_train, nb_split);
timeCls3 = toc;

% ************************************************
% results:
% timeCls = [0.0527 0.0262 0.0263];
% timeScat = [391.7722 878.4778 878.4778+17.2399]
% 
% split   1 nb train 36 accuracy 85.55 
% syn14 trans_scatt with 36 training : 85.55 += 0.00 
% split   1 nb train 36 accuracy 82.03 
% syn14 roto_trans_scatt with 36 training : 82.03 += 0.00 
% split   1 nb train 36 accuracy 82.81 
% syn14 roto_trans_scatt_log with 36 training : 82.81 += 0.00 


