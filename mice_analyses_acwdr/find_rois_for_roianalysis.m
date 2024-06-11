%% Find brain regions mice
clear, clc, close all
cd C:\Users\duodenum\Desktop\brain_stuff\mice_regularacw
mice = {'cm124', 'cm125', 'cm126', 'cm127', 'cm128'};
nmice = length(mice);
topo = cell(1, nmice);

for i = 1:nmice
    files = dir("C:\Users\duodenum\Desktop\brain_stuff\mice_wfoi\WFOM_extracted_signals\WFOM_extracted_signals\" + mice{i} + "*");
    data = load("C:\Users\duodenum\Desktop\brain_stuff\mice_wfoi\WFOM_extracted_signals\WFOM_extracted_signals\" + files(1).name);
    topo{i} = data.info.rois;
end

roidims = size(topo{1});
nroi = 92;

load("acwdr_ca.mat")

%% 
imagesc(topo{1})
colormap(jet(93))
datacursormode on
region_a = 1:8;
region_b = 9:16;
region_c = 17:23;
region_d = 24:27;
region_e = 28:34;
region_f = 35:46;

region_ap = region_a + 46;
region_bp = region_b + 46;
region_cp = region_c + 46;
region_dp = region_d + 46;
region_ep = region_e + 46;
region_fp = region_f + 46;

%% TEst
all_rois = {region_a, region_b, region_c, region_d, region_e, region_f, ...
    region_ap, region_bp, region_cp, region_dp, region_ep, region_fp};
close all

topo_test = topo{5};
for i = 1:92
    c = 1;
    for j = all_rois
        if any(i == j{1})
            topo_test(topo_test == i) = c;
        end
    c = c + 1;
    end
end

% imagesc(topo_test)

names_nice = ["Left Frontal", "Left Motor", ...
    "Left Somatosensory 1", "Left Somatosensory 2", "Left Somatosensory 3", ...
    "Left Visual", ...
    "Right Frontal", "Right Motor", ...
    "Right Somatosensory 1", "Right Somatosensory 2", "Right Somatosensory 3", ...
    "Right Visual"];
%% Visualize
cm = slanCM("Set3", 12);
i_topo = topo_test;
i_topo(i_topo == 0) = nan;
h = imagesc(i_topo);
set(h, 'AlphaData', 1 - isnan(i_topo)) % get rid of background
colormap(cm)
cb = colorbar;
axis square
axis off
set(cb, 'Ticks', 1.4:0.93:12)
set(cb, 'TickLabels', names_nice)

exportgraphics(gcf, "revision\micerois.jpg", "Resolution", 600)
%% save
save("revision\mice_rois.mat", ...
    "region_a", "region_b", "region_c", "region_d", "region_e", "region_f", ...
    "region_ap", "region_bp", "region_cp", "region_dp", "region_ep", "region_fp")