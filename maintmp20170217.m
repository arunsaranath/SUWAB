setupfile = 'C:\tetracorder4.4\mapping\mapping_cuprite\cmd.lib.setup.t4.4a5s6';
libfile = 'C:\tetracorder4.4\speclib\library06.conv\s06av95a';
lib06afile = 'C:\tetracorder4.4\speclib\library06\splib06a';
libtitlesfile = 'C:\tetracorder4.4\speclib\library06.conv\splib06b.titles';

% read library and setupfile
[recs_av95,~,rawSpcs_av95] = hyperReadSpecprYuki(libfile);
[recs,~,rawSpcs] = hyperReadSpecprYuki(lib06afile);
% read references
[refList,vals,gnotFeats,VegRatios,rawTexts] = readlibsetupall(setupfile);
refListg1 = searchby('g_id',[0,1],refList);
refListg2 = searchby('g_id',[0,2],refList);

[ spcs_av95 ] = format_s06av95a( rawSpcs_av95, recs_av95, libtitlesfile );

%% read image file
hsi_file = 'C:\tetracorder4.4\cubes\cuprite95';
hdr_file = 'C:\tetracorder4.4\cubes\cuprite95.hdr';
[hdr] = envihdrreadx(hdr_file);
[hsi] = envidataread(hsi_file,hdr);
hsi = hsi./10000;
% hsi2d = reshape(hsi,[hdr.samples*hdr.lines,hdr.bands])';

% remove noisy bands;
active_bands = [4:103 114:147 178:217];
allbands = 1:hdr.bands;
noisy_bands = setdiff(allbands,active_bands);
hsi(:,:,noisy_bands) = nan;

[wsort iwsort] = sort(hdr.wavelength,'ascend');
hsi = hsi(:,:,iwsort);

hsiRGB = double(hsi(:,:,[24,16,12]));
hsiRGB(hsiRGB>1) = 1;
hsiRGB(hsiRGB<0) = 0;

figure; ax_spc = subplot(1,1,1);
figure; imagesc(hsiRGB); set_figsize(gcf,500,600);
set(gca,'dataAspectRatio',[1 1 1]);
hdt = datacursormode(gcf);
set(hdt.CurrentDataCursor, 'Marker','+', 'MarkerFaceColor','b');
set(hdt,'UpdateFcn',{@map_cursor_default,ax_spc,hsi,wsort});

%% read classification map
res_dir = 'C:\tetracorder4.4\mapping\mapping_cuprite\group.2um';
[refListg2,summaryg2] = readResMap(refListg2,res_dir);
mapg2s = cat(3,refListg2.map);
mapg2Class = zeros(hdr.lines,hdr.samples);
mapg2Score = sum(mapg2s,3);
for i=1:length(refListg2)
    mapg2Class(mapg2s(:,:,i)>0) = i;
end

res_dir = 'C:\tetracorder4.4\mapping\mapping_cuprite\group.1um';
[refListg1,summaryg1] = readResMap(refListg1,res_dir);
mapg1s = cat(3,refListg1.map);
mapg1Class = zeros(hdr.lines,hdr.samples);
mapg1Score = sum(mapg1s,3);
for i=1:length(refListg1)
    mapg1Class(mapg1s(:,:,i)>0) = i;
end

mapScores = cat(3,mapg1Score,mapg2Score);

goodPos = sum(mapScores>50,3)>1;

%% reconstruct

s=519; l=440;
spc = squeeze(hsi(l,s,:));

% check the label
g1cls = refListg1(mapg1Class(s,l));
g2cls = refListg2(mapg2Class(s,l));

spc_g1cls = searchby('irecno',g1cls.irecno,spcs_av95);
spc_g2cls = searchby('irecno',g2cls.irecno,spcs_av95);

figure;
plot(spc_g1cls.wavelengths,spc_g1cls.reflectance);
hold on;
plot(spc_g2cls.wavelengths,spc_g2cls.reflectance,'r-');


Alib = double([spc_g1cls.reflectance spc_g2cls.reflectance]);
Alib = Alib(iwsort,:);

[Xsunsal] = sunsal(Alib(active_bands,:),spc(active_bands),'lambda',0,...
                  'positivity','yes','addone','no','verbose','yes');

ysunsal = Alib*Xsunsal;

res_sunsal = norm(spc(active_bands)-ysunsal(active_bands));

figure; plot(wsort,spc,'b-');
hold on;
plot(wsort,ysunsal,'r-');

[Xhuwacb,Z,C] = huwacb_admm2([],spc(active_bands),wsort(active_bands),'verbose','yes');
yhuwacb = Alib(active_bands,:)*Xhuwacb + C*Z;

res_huwacb = norm(spc(active_bands)-yhuwacb);

b_huwacb = nan(size(spc));
b_huwacb(active_bands) = C*Z;
spc_huwacb = nan(size(spc));
spc_huwacb(active_bands) = yhuwacb;


figure; plot(wsort,spc,'b-');
hold on;
plot(wsort,ysunsal,'r-');
plot(wsort,spc_huwacb,'-','Color',[0 0.5 0]);
plot(wsort,b_huwacb,'-','Color',[0 0.5 0]);


%%
hsi2d = reshape(hsi,[hdr.samples*hdr.lines,hdr.bands])';
[Xhuwacb,Z,C] = huwacb_admm2([],hsi2d(active_bands,:),wsort(active_bands),'verbose','yes');