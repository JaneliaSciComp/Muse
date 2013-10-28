load('localized_etc_alt.mat');
load('localized_etc_alt_jpn.mat');
alt=localized_etc_alt;
jpn=localized_etc_alt_jpn';

localized_alt=[alt.localized]';
localized_jpn=[jpn.localized]';

is_different=(localized_alt~=localized_jpn);
i_different=find(is_different);

% JPN localized 7 more than ALT on balance, but
% they differ for 105 segments of 3724 !

alt_different=alt(is_different);
jpn_different=jpn(is_different);

