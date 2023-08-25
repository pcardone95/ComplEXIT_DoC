%% Find channel location for Minji
load('chlocs_nexstim')
chn_to_find = {'Fpz','AF1','AF2','Cz','Iz'};
chn_to_find_num = [2 4 6 29 60];
sph_theta_l = [0;  11.2305; -11.2305; 90; 180;];
sph_phi_l = [-2.6980; 17.9035; 17.9035; 0; -25.86;];
eximia_addchn = chlocs(chn_to_find_num);
for chani = 1:length(chn_to_find)
    eximia_addchn(chani).labels = chn_to_find{chani};
    eximia_addchn(chani).type = [];
    eximia_addchn(chani).ref = 'average';
    eximia_addchn(chani).badchan = 1; % Dunno if necessary
    eximia_addchn(chani).urchan = 127 + chani;
    eximia_addchn(chani).sph_radius = 85;
    eximia_addchn(chani).sph_theta = sph_theta_l(chani);
    eximia_addchn(chani).sph_phi = sph_phi_l(chani);
    [eximia_addchn(chani).X, eximia_addchn(chani).Y, eximia_addchn(chani).Z] = ...
        sph2cart(deg2rad(eximia_addchn(chani).sph_theta), deg2rad(eximia_addchn(chani).sph_phi), eximia_addchn(chani).sph_radius);
end
eximia_addchn = rmfield(eximia_addchn, {'sph_theta_besa','sph_phi_besa'});
save('eximia_addchn.mat','eximia_addchn')

% Doesn't seem to be precise unfortunately at the moment.
