% This script plots spatial information that wasn't used for regrouping, as
% a sanity check to compare with existing literature on parallel fibre
% bouton information (Pichitpornchai et al., 1994; Shepherd et al., 2002)

d = get_interbouton_dist(Ain,[d1,d2],ix_axons_to_rois,Pixel_size);
[mean(d),std(d)]