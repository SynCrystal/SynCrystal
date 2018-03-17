% rsds_classif : a function to reproduce classification experiments of paper
%
%   ``Rotation, Scaling and Deformation Invariant Scattering
%   for Texture Discrimination"
%   Laurent Sifre, Stephane Mallat
%   Proc. IEEE CVPR 2013 Portland, Oregon
%

function error = rsds_classif2(db, db_name, feature_name, grid_train, nb_split)

n_class = numel(db.src.classes);
n_images = numel(db.src.files);
n_per_class = n_images / n_class; % assumes constant number of image per class

for i_split = 1:nb_split
    for i_train = 1:numel(grid_train)
        n_train = grid_train(i_train);
        train_set = [1:36,293:328,585:620,877:912]; 
        test_set = [37:292,329:584,621:876,913:1168];
        train_opt.dim = n_train;
        model = affine_train(db, train_set, train_opt);
        labels = affine_test(db, model, test_set);
        error(i_split, i_train) = classif_err(labels, test_set, db.src);
        fprintf('split %3d nb train %2d accuracy %.2f \n', ...
            i_split, n_train, 100*(1-error(i_split, i_train)));
    end
end

%% averaged performance
perf = 100*(1-mean(error));
perf_std = 100*std(error);

for i_train = 1:numel(grid_train)
    fprintf('%s %s with %2d training : %.2f += %.2f \n', ...
        db_name, feature_name, grid_train(i_train), perf(i_train), perf_std(i_train));
end

end
