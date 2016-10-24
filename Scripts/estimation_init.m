%% Initialize Estimation
% No need to edit code beyond this line

global options_

%% Distance between Bounds
% No need to edit this

% (Only used in grid search)

options_.EST.DIST = [];
FF = options_.EST.variables;
for ii = 1:length(FF)
    var_ii = char(FF(ii));
    options_.EST.DIST.(var_ii) = options_.EST.UB.(var_ii) - options_.EST.LB.(var_ii);
end

%% Add Weights to IRFs

% add weights to irfs
pvarcoirfs_clean.weight = ones(size(pvarcoirfs_clean.irf));

% rd weights
match_rd = strmatch('rd : rd', pvarcoirfs_clean.id1);
pvarcoirfs_clean(match_rd, :).weight = options_.EST.weight_rd * pvarcoirfs_clean(match_rd, :).weight;

% tfp weights
match_tfp = strmatch('rd : tfp', pvarcoirfs_clean.id1);
pvarcoirfs_clean(match_tfp, :).weight = options_.EST.weight_tfp * pvarcoirfs_clean(match_tfp, :).weight;

