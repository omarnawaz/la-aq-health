%% Documentati

% runme_v1.m

clear; clc; close all; tic;

% M. Omar Nawaz
% July 20th, 2022
% NOTES: This script generates netcdf files for input into the CAL2035
% jupyter notebook starting from the direct model outputs from the 
% GEOS-Chem adjoint. This script was written using MATALB v2020b.

%% Settings

% Path to script
base_dir = './';
% City Initials
INIT     = 'LOS';
% City Name
NAME     = 'Los Angeles';
% State city is in
STATE    = 'California';
 
%% Set-Up

% Directories
out_dir  = [base_dir,'outputs/'];
in_dir   = [base_dir,'static/'];
sens_dir = [base_dir,'inputs/',INIT,'/'];
emis_dir = [in_dir,'emissions/'];
othe_dir = [in_dir,'other/'];
ozon_dir = [in_dir,'ozone/'];
nei_dir  = [in_dir,'nei/'];
hlth_dir = [in_dir,'health/'];

% Files
pmfile     = [othe_dir,'V5GL02.HybridPM25.Global.201101-201112.nc'];
nofile     = [othe_dir,'TROP_2011_DC_US_v2.nc'];
popfile    = [othe_dir,'gpw_v4_population_count_rev11_2010_30_sec.tif'];
maskfile   = [othe_dir,'cities_13k_largerext.tif'];
lookupfile = [othe_dir,'Revised_GHS-SMOD HDC_Lookup_DAM_112421.csv'];
areafile   = [othe_dir,'GC_M2_US_05x0667.nc'];
gbdpopfile = [hlth_dir,'IHME_GBD_2019_POP_SYA_2011_Y2021M01D28.CSV'];
ctypopfile = [hlth_dir,'USCITIES_07212022_Data.csv'];
mortfile   = [hlth_dir,'IHME-GBD_2019_DATA.csv'];
asthfile   = [hlth_dir,'IHME-GBD_2019_DATA_ASTHMA.csv'];
% file containing codes for age, location, and cause
acodefile = [hlth_dir,'IHME_GBD_2019_CONTEXTS_BY_AGE_Y2020M12D18.XLSX'];
lcodefile = [hlth_dir,'IHME_GBD_2019_ALL_LOCATIONS_HIERARCHIES_Y2022M01D20.XLSX'];
ccodefile = [hlth_dir,'IHME_GBD_2019_CAUSE_HIERARCHY_Y2020M11D25.XLSX'];
files_o3  = dir([ozon_dir,'*.txt']); flen_o3 = length(files_o3);

% Assignments
IIPAR     = 91;
JJPAR     = 89;
dmax      = 396;
k_pm      = 34;
k_o3      = 62;
k_no      = 34;
ks        = [k_pm,k_o3,k_no];
di        = datenum('11/30/2010');
pollutant = {'PM','O3','NO'};
pmax      = length(pollutant);
% emission factors
SPD        = 60 * 60 * 24;
CM2M       = 10000;
MW         = [46.005, 12.0107+15.9994,12.0107];
KGCONV     = MW * 0.001 / (6.0221415E23);
efactor    = SPD * CM2M * KGCONV;
% health information
ages       = [25,30,35,40,45,50,55,60,65,70,75,80,85,90,95];
ageids     = [10,11,12,13,14,15,16,17,18,19,20,30,31,32,235];
amax       = length(ages);
% mortality rate outcome names
mortid      = [976,509,426,322,494,493];
mmax        = length(mortid);     

% Initializations
% Sensitivities
PMGDT = zeros(IIPAR,JJPAR,dmax,k_pm); PMEMS = PMGDT;
O3GDT = zeros(IIPAR,JJPAR,dmax,k_o3); O3EMS = O3GDT;
NOGDT = zeros(IIPAR,JJPAR,dmax,k_no); NOEMS = NOGDT;
% Other
inames    = cell(62,3);
O3        = zeros(91,89); CT = zeros(91,89);

%% (1) Read in adjoint sensitivities and aggregate to annual time-scale
fprintf('Begin adjoint sensitivity read and aggregation at %0.2fs\n',toc)

% loop through three pollutants
for p = 1:pmax
    A   = zeros(IIPAR,JJPAR,dmax,ks(p)); G = zeros(IIPAR,JJPAR,dmax,ks(p));
    TOT = zeros(ks(p),1); C = 0;
% loop through months    
for m = 1:12
    % get date quantities
    mon  = pad(num2str(m),2,'left','0'); date = datenum([mon,'/01/2011']);
    d0   = addtodate(date,-1,'month')-di;
    df   = addtodate(addtodate(date,+1,'month'),-1,'day')-di;
    dlen = (df-d0+1);
    % get paths to sensitivity files and skip empty months
    adj = [sens_dir,'ems.adj.' ,mon,'.',INIT,'.',pollutant{p},'.nc'];
    gdt = [sens_dir,'gctm.gdt.',mon,'.',INIT,'.',pollutant{p},'.nc'];
    if ~isfile(adj); continue; end
    I = ncinfo(adj); ivar = I.Variables; ilen = length(ivar);
    if any(TOT~=0); C = 1; end
    for i = 7:ilen        
       if C == 1; if TOT(i-6)==0; continue; end; end
       name = ivar(i).Name; inames{i-6,p} = ['IJ-GDE-S__',name(10:end)];
       tmpg = ncread(gdt,inames{i-6,p}); if sum(tmpg(:)) == 0; continue; end
       tmpa = ncread(adj,name); 
       A(:,:,d0:df,i-6) = A(:,:,d0:df,i-6) + tmpa(:,:,1:dlen) * dlen;
       G(:,:,d0:df,i-6) = G(:,:,d0:df,i-6) + tmpg(:,:,1:dlen);     
       TOT(i-6) = sum(tmpg(:));
    end
end
    if p == 1; PMEMS = A; PMGDT = G; end
    if p == 2; O3EMS = A; O3GDT = G; end
    if p == 3; NOEMS = A; NOGDT = G; end
    fprintf('Finished reading sens for %s at %0.2fs\n',pollutant{p},toc)
end

% coordinate and area data from GEOS-Chem
latgc = ncread(adj,'LAT'); 
longc = ncread(adj,'LON');
GCM2  = ncread(areafile,'DXYP__DXYP');

%% (2) Get cost-function data
fprintf('Beginning cost-function calculation %0.2fs\n',toc)

% Read in mask data
%-------------------------------------------------------------------------%
[A,~]     = readgeoraster(maskfile,'CoordinateSystemType','planar');
T         = readtable(lookupfile,'VariableNamingRule','preserve'); 
ind       = T{find(strcmp(T{:,2},NAME)),1};
A(A~=ind) = 0;
% coordinates
latf = linspace(75,-56,15720); lonf = linspace(-179.5000,179.0000,43020);
% get limits
xmin = find(sum(A,1)>0,1,'first'); xmax = find(sum(A,1)>0,1,'last');
ymin = find(sum(A,2)>0,1,'first'); ymax = find(sum(A,2)>0,1,'last');
% limit to just the city boundaries
Z  = A(ymin:ymax,xmin:xmax); la = latf(ymin:ymax); lo = lonf(xmin:xmax);
fprintf('Finished reading in mask data %0.2fs\n',toc)
%-------------------------------------------------------------------------%

% Read in O3 data
%-------------------------------------------------------------------------%
for f = 1:flen_o3
    temp = load([ozon_dir,files_o3(f).name])';
    O3(temp>0) = O3(temp>0) + temp(temp>0); CT(temp>0) = CT(temp>0) + 1;
end
O3 = O3 ./ CT;
fprintf('Finished reading in o3 data %0.2fs\n',toc)
%-------------------------------------------------------------------------%

% Read in PM Data
%-------------------------------------------------------------------------%
PM       = ncread(pmfile,'GWRPM25'); 
latpm    = ncread(pmfile,'lat'); lonpm = ncread(pmfile,'lon');
fprintf('Finished reading in pm2.5 data %0.2fs\n',toc)
%-------------------------------------------------------------------------%

% Read in Pop data
%-------------------------------------------------------------------------%
[POP, R] = readgeoraster(popfile); POP = flipud(POP);
latpo    = linspace(-90,89.999999,21600); lonpo = linspace(-180,180,43200);
fprintf('Finished reading in population data %0.2fs\n',toc)
%-------------------------------------------------------------------------%

% Read in NO2 data
%-------------------------------------------------------------------------%
NO       = ncread(nofile,'DS'); 
latno    = ncread(nofile,'lat'); lonno = ncread(nofile,'lon');
fprintf('Finished reading in no2 data %0.2fs\n',toc)
%-------------------------------------------------------------------------%

% Calculate city mask at 0.5° x 0.667° resolution
%-------------------------------------------------------------------------%
CITY_MASK = zeros(91,89);
for i = 1:length(lo)
    for j = 1:length(la)
        [~,x] = min(abs(longc-lo(i))); [~,y] = min(abs(latgc-la(j)));        
        if Z(j,i) ~= 0; CITY_MASK(x,y) = CITY_MASK(x,y) + 1/4840; end
    end
end
%-------------------------------------------------------------------------%

% Calculate state mask at at 0.5° x 0.667° resolution
%-------------------------------------------------------------------------%
% read in state shape file
S = shaperead([othe_dir,'s_22mr22.shp']);
% loop through states until matching state
for s = 1:59; if strcmp(S(s).NAME,STATE); break; end; end
% get polygon vertices and create meshgrid of geos chem coordinates
XS = S(s).X; YS = S(s).Y; [X,Y] = meshgrid(longc,latgc);
% coarse mask
in = double(inpolygon(X,Y,XS,YS));
% get fractional grid cells values of coarse cells
for I = 2:IIPAR-1
for J = 2:JJPAR-1
    C = 0; dJ = (J-1):(J+1); dI = (I-1):(I+1); n = 2; err = 1;
    if (in(J,I) == 0) && any(in(dJ,dI)==1,'all'); C = 1; end
    if (in(J,I) == 1) && any(in(dJ,dI)==0,'all'); C = 1; end
    if C == 0; continue; end
    while err > 0.05
        ln      = linspace((longc(I)-1/3),(longc(I)+1/3),n);
        lt      = linspace((latgc(J)-1/4),(latgc(J)+1/4),n);
        [X,Y]   = meshgrid(ln,lt);
        tmp     = sum(double(inpolygon(X,Y,XS,YS)),'all');
        err     = abs(in(J,I) - tmp / n^2);
        in(J,I) = tmp / n^2;
        n       = n + 2;        
    end
end    
end
STATE_MASK = in';
%-------------------------------------------------------------------------%


% Save out mask files
%-------------------------------------------------------------------------%
fname = [out_dir,'CITY_05x06.csv']; if isfile(fname); delete(fname); end
csvwrite(fname,CITY_MASK');
 
fname = [out_dir,'STATE_05x06.csv']; if isfile(fname); delete(fname); end
csvwrite(fname,STATE_MASK');
fprintf('Saved coarse mask files %0.2fs\n',toc) 
%-------------------------------------------------------------------------%

%% (3) Calculate cost-function values

% PM 2.5
NUM = 0; DEN = 0;
for i = 1:length(lo)
    for j = 1:length(la)
        [~,x] = min(abs(lonpm-lo(i)));
        [~,y] = min(abs(latpm-la(j)));
        [~,xp] = min(abs(lonpo-lo(i)));
        [~,yp] = min(abs(latpo-la(j)));        
        if isnan(PM(x,y)); continue; end
        if POP(yp,xp) <=0; continue; end
        NUM = NUM + double(PM(x,y)) * double(Z(j,i)) * double(POP(yp,xp));
        DEN = DEN + double(Z(j,i))  * double(POP(yp,xp));
    end
end
CFPM = NUM / DEN;
fprintf('Finished calculating PM2.5 cost-function %0.2fs\n',toc)

% O3
NUM = 0; DEN = 0;
for i = 1:length(lo)
    for j = 1:length(la)
        [~,x] = min(abs(longc-lo(i)));
        [~,y] = min(abs(latgc-la(j)));
        [~,xp] = min(abs(lonpo-lo(i)));
        [~,yp] = min(abs(latpo-la(j)));        
        if isnan(O3(x,y)); continue; end
        if POP(yp,xp) <=0; continue; end
        NUM = NUM + double(O3(x,y)) * double(Z(j,i)) * double(POP(yp,xp));
        DEN = DEN + double(Z(j,i))  * double(POP(yp,xp));
    end
end
CFO3 = NUM / DEN;
fprintf('Finished calculating O3 cost-function %0.2fs\n',toc)

% NO2
NUM = 0; DEN = 0;
for i = 1:length(lo)
    for j = 1:length(la)
        [~,x] = min(abs(lonno-lo(i)));
        [~,y] = min(abs(latno-la(j)));
        [~,xp] = min(abs(lonpo-lo(i)));
        [~,yp] = min(abs(latpo-la(j)));        
        if isnan(NO(x,y)); continue; end
        if POP(yp,xp) <=0; continue; end
        NUM = NUM + double(NO(x,y)) * double(Z(j,i)) * double(POP(yp,xp));
        DEN = DEN + double(Z(j,i))  * double(POP(yp,xp));
    end
end
CFNO = NUM / DEN;
fprintf('Finished calculating NO2 cost-function %0.2fs\n',toc)

%% (4) Emission scaling for gctm.gdt files

for d = 1:dmax
    % establish date
    date    = char(datetime(2010,11,30,'Format','yyyyMMdd') + days(d)); 
    date    = strrep(date,'2010','2011'); 
    vocfile = [emis_dir,'ctm.VOC.',date,'.nc'];
    noxfile = [emis_dir,'ctm.NOx.',date,'.nc'];
    % NOx scaling
    ems   = sum(ncread(noxfile,'NOX-AN-S__NOx'),3) * efactor(1) .* GCM2 ; 
    PMEMS(:,:,d,25) = PMGDT(:,:,d,25) ./ ems;
    O3EMS(:,:,d,25) = O3GDT(:,:,d,25) ./ ems;
    NOEMS(:,:,d,25) = NOGDT(:,:,d,25) ./ ems;    
    %VOC scaling for O3
    for k = 1:k_o3
        kname = erase(erase(inames{k,2},'IJ-GDE-S__'),'_an');
        try
            temp  = ncread(vocfile,['ANTHSRCE__',kname]);
        catch
            continue
        end
        if strcmp(kname,'CO')
            ems = sum(temp,3) * efactor(2) .* GCM2;
        else
            ems = sum(temp,3) * efactor(3) .* GCM2;            
        end        
        cf = mean(mean(ems(ems>0))) / (1E1);           
        tmp = O3GDT(:,:,d,k) ./ ems; tmp(ems < cf) = 0; O3EMS(:,:,d,k) = tmp;         
    end
end

% remove infinities and NaNs from emission division
PMEMS( PMEMS ==  Inf  ) = 0;
PMEMS( PMEMS == -Inf  ) = 0;
PMEMS( isnan( PMEMS ) ) = 0;
O3EMS( O3EMS ==  Inf  ) = 0;
O3EMS( O3EMS == -Inf  ) = 0;
O3EMS( isnan( O3EMS ) ) = 0;
NOEMS( NOEMS ==  Inf  ) = 0;
NOEMS( NOEMS == -Inf  ) = 0;
NOEMS( isnan( NOEMS ) ) = 0;
fprintf('Finished emission scaling of sensitivities %0.2fs\n',toc)

%% (5) Apply cost-function scaling

% PM2.5
sf    = CFPM / sum(PMGDT,'all'); 
PMEMS = PMEMS * sf;
fprintf('Scaling factor of %0.2f applied for pm2.5 at %0.2fs\n',sf,toc)
% O3
sf    = CFO3 / sum(O3GDT,'all'); 
O3EMS = O3EMS * sf;
fprintf('Scaling factor of %0.2f applied for o3 at %0.2fs\n',sf,toc)
% NO2
sf    = CFNO / sum(NOGDT,'all'); 
NOEMS = NOEMS * sf;
fprintf('Scaling factor of %0.2f applied for no2 at %0.2fs\n',sf,toc)

%% (6) Calculate contributions

% Read in NEI emissions
files   = dir([nei_dir,'nei_emis*.nc']);
flen = length(files);
E     = zeros(91,89,396,k_o3);
for f = 1:flen
    pathname = [files(f).folder,'/',files(f).name];
    for i = 1:k_o3
        if contains(inames{i,2},'_an')
            if contains(inames{i,2},'PO'); continue; end
            iname = erase(inames{i,2}(11:end),{'_an1','_an','PI'});
            if contains(inames{i,2},'an2'); iname = 'SO22'; end            
            if strcmp(iname,'NOX'); iname = 'NOx'; end
            if strcmp(iname,'ISOP'); continue; end
            species{i} = iname;
            E(:,:,:,i) = E(:,:,:,i) + ncread(pathname,iname);
        end
    end 
end

% Calculate contributions
for i = 1:62
    for p = 1:3
        k = find(strcmp(inames(i,2),inames(:,p)));
        if isempty(k); continue; end
        species_f{k,p} = species{i};
        if p == 1; dJ_PM(:,:,:,k) = E(:,:,:,i) .* PMEMS(:,:,:,k); end
        if p == 2; dJ_O3(:,:,:,k) = E(:,:,:,i) .* O3EMS(:,:,:,k); end
        if p == 3; dJ_NO(:,:,:,k) = E(:,:,:,i) .* NOEMS(:,:,:,k); end
    end
end

fprintf('Finished exposure contribution calculation %0.2fs\n',toc)

%% (7) Preprocess health data

% Load population data, mortality rates, and GBD codes
CPOP    = readtable(ctypopfile); 
cind    = find(strcmp(CPOP{:,3},NAME));
CITYPOP = CPOP{cind,7}; 
POP_PED = CPOP{cind,23};
GBDPOP  = readtable(gbdpopfile); 
GBDPOP  = GBDPOP(find(strcmp('both',GBDPOP{:,4})),:);
sind    = find(strcmp(GBDPOP{:,2},STATE));
GBDPOP  = GBDPOP(sind,:);
MORT    = readtable(mortfile); ASTH = readtable(asthfile);
A       = readtable(acodefile); 
L       = readtable(lcodefile); 
C       = readtable(ccodefile);

% Process population data
POP = zeros(amax,1);
LB  = zeros(amax,1);
UB  = zeros(amax,1);
% loop through age groups and calculate pop in age brackets
for a = 1:amax       
   aind = ((GBDPOP{:,6} >= ages(a)) & (GBDPOP{:,6} < ages(a)+5));
   POP(a)=sum(GBDPOP{aind,12});
   LB(a)=sum(GBDPOP{aind,14})/sum(GBDPOP{aind,12}); 
   UB(a)=sum(GBDPOP{aind,13})/sum(GBDPOP{aind,12}); 
   if a == amax
       temp = GBDPOP{:,12};
       POP(a) = temp(end); 
       temp = GBDPOP{:,14};       
       LB(a) = temp(end) / POP(a);
       temp = GBDPOP{:,13};       
       UB(a) = temp(end) / POP(a);       
   end
end
% get all ages population
TOT     = sum(GBDPOP{:,12});   
TOT_LB  = sum(GBDPOP{:,14});   
TOT_UB  = sum(GBDPOP{:,13});   
% calculate bounds of PED pop
POP_PED_LB = POP_PED * (TOT_LB -sum(POP.*LB)) / (TOT-sum(POP));
POP_PED_UB = POP_PED * (TOT_UB -sum(POP.*UB)) / (TOT-sum(POP));
% calculate final city pop
POP  = CITYPOP * (POP / TOT);
POP_LB = LB .* POP;
POP_UB = UB .* POP;

% Process mortality rates
lid       = L{find(strcmp(L{:,3},STATE)),2};
rate      = zeros(amax,mmax);    
rate_lb   = zeros(amax,mmax);    
rate_ub   = zeros(amax,mmax);    
rate_asth = ASTH{find(lid == ASTH{:,2}),8};
rate_asth_lb = ASTH{find(lid == ASTH{:,2}),10};
rate_asth_ub = ASTH{find(lid == ASTH{:,2}),9};
% loop through mortality rates
for m = 1:mmax
   mid   = mortid(m);
% loop through ages
for a = 1:amax
    aid = ageids(a);
    inds = find((MORT{:,4} == aid) & (MORT{:,2} == lid) & (MORT{:,5} == mid));
    rate(a,m) = MORT{inds,8};
    rate_lb(a,m) = MORT{inds,10};
    rate_ub(a,m) = MORT{inds,9};
end   
end

fprintf('Finished processing health data %0.2fs\n',toc)

%% (8) Health impact contribution
 
% PM2.5
%-------------------------------------------------------------------------%
% load city-independent health data
load([hlth_dir,'rr.mat']); 
load([hlth_dir,'rr_lb.mat']); 
load([hlth_dir,'rr_ub.mat']); 
load([hlth_dir,'exp.mat'])
% initialize outputs
dJPMD    = zeros(IIPAR,JJPAR,dmax,k_pm);
dJPMD_lb = zeros(IIPAR,JJPAR,dmax,k_pm);
dJPMD_ub = zeros(IIPAR,JJPAR,dmax,k_pm);
AVAIL  = squeeze(sum(dJ_PM,1:3));
% calculate health impacts
for k = 1:k_pm
    if AVAIL(k) == 0; continue; end
    [base_pm,dJPMD(:,:,:,k)]       = pm_health_calc(POP,rate,EXP, ...
        RR,CFPM,dJ_PM(:,:,:,k));
    [base_pm_lb,dJPMD_lb(:,:,:,k)] = pm_health_calc(POP_LB,rate_lb,EXP, ...
        RR_LB,CFPM,dJ_PM(:,:,:,k));
    [base_pm_ub,dJPMD_ub(:,:,:,k)] = pm_health_calc(POP_UB,rate_ub,EXP, ...
        RR_UB,CFPM,dJ_PM(:,:,:,k));
end
fprintf('Finished PM2.5 health contribution calculation %0.2fs\n',toc)
%-------------------------------------------------------------------------%

% O3
%-------------------------------------------------------------------------%
% initialize outputs
dJO3D    = zeros(IIPAR,JJPAR,dmax,k_o3);
dJO3D_lb = zeros(IIPAR,JJPAR,dmax,k_o3);
dJO3D_ub = zeros(IIPAR,JJPAR,dmax,k_o3);
AVAIL    = squeeze(sum(dJ_O3,1:3));
beta     = log(1.06) / 10;
beta_lb  = log(1.03) / 10;
beta_ub  = log(1.10) / 10;
tmrel    = 32.4;
tmrel_lb = 35.7;
tmrel_ub = 29.1;
% calculate health impacts
for k = 1:k_pm
    if AVAIL(k) == 0; continue; end         
    [base_o3,dJO3D(:,:,:,k)] = o3_health_calc(POP, rate(:,2), ...
        CFO3, dJ_O3(:,:,:,k), beta, tmrel);  
    [base_o3_lb,dJO3D_lb(:,:,:,k)] = o3_health_calc(POP_LB, rate_lb(:,2), ...
        CFO3, dJ_O3(:,:,:,k), beta_lb, tmrel_lb);  
    [base_o3_ub,dJO3D_ub(:,:,:,k)] = o3_health_calc(POP_UB, rate_ub(:,2), ...
        CFO3, dJ_O3(:,:,:,k), beta_ub, tmrel_ub);  
end
fprintf('Finished O3 health contribution calculation %0.2fs\n',toc)
%-------------------------------------------------------------------------%

% NO2
%-------------------------------------------------------------------------%

% initialize outputs
dJNOD    = zeros(IIPAR,JJPAR,dmax,k_no);   
dJNOD_lb = zeros(IIPAR,JJPAR,dmax,k_no);   
dJNOD_ub = zeros(IIPAR,JJPAR,dmax,k_no);   
AVAIL    = squeeze(sum(dJ_NO,1:3));
beta     = log(1.26) / 10;
beta_lb  = log(1.10) / 10;
beta_ub  = log(1.37) / 10;
tmrel    = 2;
tmrel_lb = 5;
tmrel_ub = 0;
for k = 1:k_no
    if AVAIL(k) == 0; continue; end     
    [base_no,dJNOD(:,:,:,k)] = no_health_calc(POP_PED, rate_asth, ...
        CFNO,dJ_NO(:,:,:,k),beta,tmrel);  
    [base_no_lb,dJNOD_lb(:,:,:,k)] = no_health_calc(POP_PED_LB, rate_asth_lb, ...
            CFNO,dJ_NO(:,:,:,k),beta_lb,tmrel_lb);      
    [base_no_ub,dJNOD_ub(:,:,:,k)] = no_health_calc(POP_PED_UB, rate_asth_ub, ...
            CFNO,dJ_NO(:,:,:,k),beta_ub,tmrel_ub);             
end

fprintf('Finished NO2 health contribution calculation %0.2fs\n',toc)
%-------------------------------------------------------------------------%

% Calculate health scaling factors
PMHF = dJPMD ./ dJ_PM; PMHF_LB = dJPMD_lb ./ dJ_PM; PMHF_UB = dJPMD_ub ./ dJ_PM;
O3HF = dJO3D ./ dJ_O3; O3HF_LB = dJO3D_lb ./ dJ_O3; O3HF_UB = dJO3D_ub ./ dJ_O3;
NOHF = dJNOD ./ dJ_NO; NOHF_LB = dJNOD_lb ./ dJ_NO; NOHF_UB = dJNOD_ub ./ dJ_NO;
% Remove NaN values
PMHF(isnan(PMHF)) = 0; PMHF_LB(isnan(PMHF_LB)) = 0; PMHF_UB(isnan(PMHF_UB)) = 0;
O3HF(isnan(O3HF)) = 0; O3HF_LB(isnan(O3HF_LB)) = 0; O3HF_UB(isnan(O3HF_UB)) = 0;
NOHF(isnan(NOHF)) = 0; NOHF_LB(isnan(NOHF_LB)) = 0; NOHF_UB(isnan(NOHF_UB)) = 0;


%% (9) Write netCDF outputs for sensitivities

%-----------------------
% setup CAL2035 info
%-----------------------

% calculate date
date = strrep(strrep(strrep(char(datetime('now','Format',"dd/MM/uuuu HH:mm:ss")),'/',''),':',''),' ','_');

ks = {k_pm,k_o3,k_no};
% delte previous files in folder
delete([out_dir,'*.nc'])

% loop through pollutants
for p = 1:3

% output data file
fname_out = [out_dir,INIT,'_sens_',pollutant{p},'_',date,'.nc'];

% create file
ncid = netcdf.create(fname_out,'netcdf4'); 

% define dimensions
latid    = netcdf.defDim(ncid,'lat',JJPAR);
lonid    = netcdf.defDim(ncid,'lon',IIPAR);  
dayid    = netcdf.defDim(ncid,'time',dmax);

% define 1d variables
lonvarid = netcdf.defVar(ncid,'lon','NC_DOUBLE',lonid);
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');
netcdf.putAtt(ncid,lonvarid,'long_name','longitude')
netcdf.putAtt(ncid,lonvarid,'comment','centre of grid cell')

latvarid = netcdf.defVar(ncid,'lat','NC_DOUBLE',latid);
netcdf.putAtt(ncid,latvarid,'units','degrees_north');
netcdf.putAtt(ncid,latvarid,'long_name','latitude')
netcdf.putAtt(ncid,latvarid,'comment','centre of grid cell')   

dayvarid = netcdf.defVar(ncid,'time','NC_DOUBLE',dayid);
netcdf.putAtt(ncid,dayvarid,'units','days since 1981-01-01 00:00:00');
netcdf.putAtt(ncid,dayvarid,'long_name','reference time of sst field')
netcdf.putAtt(ncid,dayvarid,'comment','')  
netcdf.putAtt(ncid,dayvarid,'axis','T')

% loop through species
for k = 1:ks{p}
    spec = species_f{k,p}; if isempty(spec); continue; end

    % 4d data variable info
    NEI11_input.n  = spec;
    var4d_max      = length(NEI11_input);

    % 4d long name
    NEI11_input_name(1).n = ['daily sensitivity of ',pollutant{p},' to ',spec];

    % 4d units
    NEI11_input_unit(1 ).n = ['(ppbv) / (kg / m2 / s)'];

    % define 4d data variables
    var4d_k_id = netcdf.defVar(ncid,NEI11_input.n,'NC_DOUBLE',[lonid latid dayid]);
    netcdf.putAtt(ncid,var4d_k_id,'long_name',NEI11_input_name.n);
    netcdf.putAtt(ncid,var4d_k_id,'units',NEI11_input_unit.n); 

    % deflate variable
    netcdf.defVarDeflate(ncid,var4d_k_id, true, true, 7)   

    %  done defining variables
    netcdf.endDef(ncid)

    % input matlab arrays into file
    if p == 1; netcdf.putVar(ncid,var4d_k_id,PMEMS(:,:,:,k)); end
    if p == 2; netcdf.putVar(ncid,var4d_k_id,O3EMS(:,:,:,k)); end
    if p == 3; netcdf.putVar(ncid,var4d_k_id,NOEMS(:,:,:,k)); end

end

% add lat, lon to file
netcdf.putVar(ncid,lonvarid,longc);
netcdf.putVar(ncid,latvarid,latgc);
netcdf.putVar(ncid,dayvarid,(1:dmax) + 10925)

% file description
descr = ['This file contains sensitivity data calculated from the GEOS-Chem', ...
         ' adjoint. We calculate the sensitivity of ',spec,' with respect to ', ...
         ' precursor species at the 0.5° x 0.667° resolution', ...
         ' from December 1st 2010 to December 31st 2011.'];

% write global data
fileattrib(fname_out,'+w');
ncwriteatt(fname_out,'/','file creation time:',datestr(now));
ncwriteatt(fname_out,'/','description:',descr);

% close this file
netcdf.close(ncid);   

end

%% (10) Write netCDF outputs for health fractions

%-----------------------
% setup CAL2035 info
%-----------------------

% calculate date
date = strrep(strrep(strrep(char(datetime('now','Format',"dd/MM/uuuu HH:mm:ss")),'/',''),':',''),' ','_');

ks = {k_pm,k_o3,k_no};

% loop through pollutants
for p = 1:3

% output data file
fname_out = [out_dir,INIT,'_health_',pollutant{p},'_',date,'.nc'];

% create file
ncid = netcdf.create(fname_out,'netcdf4'); 

% define dimensions
latid    = netcdf.defDim(ncid,'lat',JJPAR);
lonid    = netcdf.defDim(ncid,'lon',IIPAR);  
dayid    = netcdf.defDim(ncid,'time',dmax);

% define 1d variables
lonvarid = netcdf.defVar(ncid,'lon','NC_DOUBLE',lonid);
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');
netcdf.putAtt(ncid,lonvarid,'long_name','longitude')
netcdf.putAtt(ncid,lonvarid,'comment','centre of grid cell')

latvarid = netcdf.defVar(ncid,'lat','NC_DOUBLE',latid);
netcdf.putAtt(ncid,latvarid,'units','degrees_north');
netcdf.putAtt(ncid,latvarid,'long_name','latitude')
netcdf.putAtt(ncid,latvarid,'comment','centre of grid cell')   

dayvarid = netcdf.defVar(ncid,'time','NC_DOUBLE',dayid);
netcdf.putAtt(ncid,dayvarid,'units','days since 1981-01-01 00:00:00');
netcdf.putAtt(ncid,dayvarid,'long_name','reference time of sst field')
netcdf.putAtt(ncid,dayvarid,'comment','')  
netcdf.putAtt(ncid,dayvarid,'axis','T')

% loop through species
for k = 1:ks{p}
    spec = species_f{k,p}; if isempty(spec); continue; end

    % 4d data variable info
    NEI11_input.n  = spec;
    var4d_max      = length(NEI11_input);

    % 4d long name
    NEI11_input_name(1).n = ['daily sensitivity of ',pollutant{p},' to ',spec];

    % 4d units
    NEI11_input_unit(1 ).n = ['(ppbv) / (kg / m2 / s)'];

    % define 4d data variables
    var4d_k_id = netcdf.defVar(ncid,NEI11_input.n,'NC_DOUBLE',[lonid latid dayid]);
    netcdf.putAtt(ncid,var4d_k_id,'long_name',NEI11_input_name.n);
    netcdf.putAtt(ncid,var4d_k_id,'units',NEI11_input_unit.n); 

    % deflate variable
    netcdf.defVarDeflate(ncid,var4d_k_id, true, true, 7)   

    %  done defining variables
    netcdf.endDef(ncid)

    % input matlab arrays into file
    if p == 1; netcdf.putVar(ncid,var4d_k_id,PMHF(:,:,:,k)); end
    if p == 2; netcdf.putVar(ncid,var4d_k_id,O3HF(:,:,:,k)); end
    if p == 3; netcdf.putVar(ncid,var4d_k_id,NOHF(:,:,:,k)); end

end

% add lat, lon to file
netcdf.putVar(ncid,lonvarid,longc);
netcdf.putVar(ncid,latvarid,latgc);
netcdf.putVar(ncid,dayvarid,(1:dmax) + 10925)

% file description
descr = ['This file contains sensitivity data calculated from the GEOS-Chem', ...
         ' adjoint. We calculate the sensitivity of ',spec,' with respect to ', ...
         ' precursor species at the 0.5° x 0.667° resolution', ...
         ' from December 1st 2010 to December 31st 2011.'];

% write global data
fileattrib(fname_out,'+w');
ncwriteatt(fname_out,'/','file creation time:',datestr(now));
ncwriteatt(fname_out,'/','description:',descr);

% close this file
netcdf.close(ncid);   

end

%% (11) Write netCDF outputs for health fractions lower bounds

%-----------------------
% setup CAL2035 info
%-----------------------

% calculate date
date = strrep(strrep(strrep(char(datetime('now','Format',"dd/MM/uuuu HH:mm:ss")),'/',''),':',''),' ','_');

ks = {k_pm,k_o3,k_no};

% loop through pollutants
for p = 1:3

% output data file
fname_out = [out_dir,INIT,'_health_lb_',pollutant{p},'_',date,'.nc'];

% create file
ncid = netcdf.create(fname_out,'netcdf4'); 

% define dimensions
latid    = netcdf.defDim(ncid,'lat',JJPAR);
lonid    = netcdf.defDim(ncid,'lon',IIPAR);  
dayid    = netcdf.defDim(ncid,'time',dmax);

% define 1d variables
lonvarid = netcdf.defVar(ncid,'lon','NC_DOUBLE',lonid);
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');
netcdf.putAtt(ncid,lonvarid,'long_name','longitude')
netcdf.putAtt(ncid,lonvarid,'comment','centre of grid cell')

latvarid = netcdf.defVar(ncid,'lat','NC_DOUBLE',latid);
netcdf.putAtt(ncid,latvarid,'units','degrees_north');
netcdf.putAtt(ncid,latvarid,'long_name','latitude')
netcdf.putAtt(ncid,latvarid,'comment','centre of grid cell')   

dayvarid = netcdf.defVar(ncid,'time','NC_DOUBLE',dayid);
netcdf.putAtt(ncid,dayvarid,'units','days since 1981-01-01 00:00:00');
netcdf.putAtt(ncid,dayvarid,'long_name','reference time of sst field')
netcdf.putAtt(ncid,dayvarid,'comment','')  
netcdf.putAtt(ncid,dayvarid,'axis','T')

% loop through species
for k = 1:ks{p}
    spec = species_f{k,p}; if isempty(spec); continue; end

    % 4d data variable info
    NEI11_input.n  = spec;
    var4d_max      = length(NEI11_input);

    % 4d long name
    NEI11_input_name(1).n = ['daily sensitivity of ',pollutant{p},' to ',spec];

    % 4d units
    NEI11_input_unit(1 ).n = ['(ppbv) / (kg / m2 / s)'];

    % define 4d data variables
    var4d_k_id = netcdf.defVar(ncid,NEI11_input.n,'NC_DOUBLE',[lonid latid dayid]);
    netcdf.putAtt(ncid,var4d_k_id,'long_name',NEI11_input_name.n);
    netcdf.putAtt(ncid,var4d_k_id,'units',NEI11_input_unit.n); 

    % deflate variable
    netcdf.defVarDeflate(ncid,var4d_k_id, true, true, 7)   

    %  done defining variables
    netcdf.endDef(ncid)

    % input matlab arrays into file
    if p == 1; netcdf.putVar(ncid,var4d_k_id,PMHF_LB(:,:,:,k)); end
    if p == 2; netcdf.putVar(ncid,var4d_k_id,O3HF_LB(:,:,:,k)); end
    if p == 3; netcdf.putVar(ncid,var4d_k_id,NOHF_LB(:,:,:,k)); end

end

% add lat, lon to file
netcdf.putVar(ncid,lonvarid,longc);
netcdf.putVar(ncid,latvarid,latgc);
netcdf.putVar(ncid,dayvarid,(1:dmax) + 10925)

% file description
descr = ['This file contains sensitivity data calculated from the GEOS-Chem', ...
         ' adjoint. We calculate the sensitivity of ',spec,' with respect to ', ...
         ' precursor species at the 0.5° x 0.667° resolution', ...
         ' from December 1st 2010 to December 31st 2011.'];

% write global data
fileattrib(fname_out,'+w');
ncwriteatt(fname_out,'/','file creation time:',datestr(now));
ncwriteatt(fname_out,'/','description:',descr);

% close this file
netcdf.close(ncid);   

end

%% (12) Write netCDF outputs for health fractions upper bounds

%-----------------------
% setup CAL2035 info
%-----------------------

% calculate date
date = strrep(strrep(strrep(char(datetime('now','Format',"dd/MM/uuuu HH:mm:ss")),'/',''),':',''),' ','_');

ks = {k_pm,k_o3,k_no};

% loop through pollutants
for p = 1:3

% output data file
fname_out = [out_dir,INIT,'_health_ub_',pollutant{p},'_',date,'.nc'];

% create file
ncid = netcdf.create(fname_out,'netcdf4'); 

% define dimensions
latid    = netcdf.defDim(ncid,'lat',JJPAR);
lonid    = netcdf.defDim(ncid,'lon',IIPAR);  
dayid    = netcdf.defDim(ncid,'time',dmax);

% define 1d variables
lonvarid = netcdf.defVar(ncid,'lon','NC_DOUBLE',lonid);
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');
netcdf.putAtt(ncid,lonvarid,'long_name','longitude')
netcdf.putAtt(ncid,lonvarid,'comment','centre of grid cell')

latvarid = netcdf.defVar(ncid,'lat','NC_DOUBLE',latid);
netcdf.putAtt(ncid,latvarid,'units','degrees_north');
netcdf.putAtt(ncid,latvarid,'long_name','latitude')
netcdf.putAtt(ncid,latvarid,'comment','centre of grid cell')   

dayvarid = netcdf.defVar(ncid,'time','NC_DOUBLE',dayid);
netcdf.putAtt(ncid,dayvarid,'units','days since 1981-01-01 00:00:00');
netcdf.putAtt(ncid,dayvarid,'long_name','reference time of sst field')
netcdf.putAtt(ncid,dayvarid,'comment','')  
netcdf.putAtt(ncid,dayvarid,'axis','T')

% loop through species
for k = 1:ks{p}
    spec = species_f{k,p}; if isempty(spec); continue; end

    % 4d data variable info
    NEI11_input.n  = spec;
    var4d_max      = length(NEI11_input);

    % 4d long name
    NEI11_input_name(1).n = ['daily sensitivity of ',pollutant{p},' to ',spec];

    % 4d units
    NEI11_input_unit(1 ).n = ['(ppbv) / (kg / m2 / s)'];

    % define 4d data variables
    var4d_k_id = netcdf.defVar(ncid,NEI11_input.n,'NC_DOUBLE',[lonid latid dayid]);
    netcdf.putAtt(ncid,var4d_k_id,'long_name',NEI11_input_name.n);
    netcdf.putAtt(ncid,var4d_k_id,'units',NEI11_input_unit.n); 

    % deflate variable
    netcdf.defVarDeflate(ncid,var4d_k_id, true, true, 7)   

    %  done defining variables
    netcdf.endDef(ncid)

    % input matlab arrays into file
    if p == 1; netcdf.putVar(ncid,var4d_k_id,PMHF_UB(:,:,:,k)); end
    if p == 2; netcdf.putVar(ncid,var4d_k_id,O3HF_UB(:,:,:,k)); end
    if p == 3; netcdf.putVar(ncid,var4d_k_id,NOHF_UB(:,:,:,k)); end

end

% add lat, lon to file
netcdf.putVar(ncid,lonvarid,longc);
netcdf.putVar(ncid,latvarid,latgc);
netcdf.putVar(ncid,dayvarid,(1:dmax) + 10925)

% file description
descr = ['This file contains sensitivity data calculated from the GEOS-Chem', ...
         ' adjoint. We calculate the sensitivity of ',spec,' with respect to ', ...
         ' precursor species at the 0.5° x 0.667° resolution', ...
         ' from December 1st 2010 to December 31st 2011.'];

% write global data
fileattrib(fname_out,'+w');
ncwriteatt(fname_out,'/','file creation time:',datestr(now));
ncwriteatt(fname_out,'/','description:',descr);

% close this file
netcdf.close(ncid);   

end


%% (13) Functions

function [base,D] = pm_health_calc(pop,rate,exps,rr,cfn,dJ)

% Function "pm_health_calc" takes inputs of population data (pop),
% mortality rates (rate), a base cost-function value (cfn), and an adjoint 
% contribution array (dJ) and calculates the total and source apportionment of
% heatlh impacts.

% Variables:
%-------------------------------------------------------------------------%
% pop  - m x 1 vector of age-stratified population data
% rate - m x n array of outcome and age seperated mortality rates
% exps - m x n cell array of exposure values from look-up tables
% rr   - m x n cell array of relative risks values from look-up tables
% cfn  - scalar cost-function value of pm2.5 exposure
% dJ   - 2D, 3D, or, 4D array of sources that contribute to the cost-function 
%        exposure

% Outputs
%-------------------------------------------------------------------------%
% base - total pm2.5 exposure deaths from cost-function
% D    - an array of source apportionment deaths

% Set-up
%-------------------------------------------------------------------------%

% Assignments
amax = length(pop);    % number of age groups
omax = size(rate,2); % number of health outcomes
% get size of source apportionment array
numVars = length(size(dJ)); sz = zeros(numVars,1);
for n = 1:numVars
   sz(n) = size(dJ,n); 
end

% Initializations
D    = zeros(sz');
base = zeros(amax,omax);

% Calculate premature deaths
%-------------------------------------------------------------------------%

% loop through age groups and 
for a = 1:amax
for o = 1:omax
     
        % extract relative risk value for age and outcome
        RR = table2array(rr{a,o}); EXP = table2array(exps{a,o});       
        % get age population
        POP = pop(a);
        % get mortality rate for outcome and age group
        y0 = rate(a,o) / 1E5;
        % get two closest look up points and exposure levels
        ind_1 = find(cfn >= EXP,1,'last'); 
        ind_2 = find(cfn <  EXP,1,'first'); 
        z0    = EXP(ind_1); zf    = EXP(ind_2);             
        % calculate fraction between points for linear extrapolation
        frac  = (cfn - z0)/(zf-z0);            
        % get lower (ind_1) and upper (ind_2) relative risks
        RR_1 = RR(ind_1); RR_2 = RR(ind_2);
        
        % calculate cost-function premature deaths for both risks
        base_d1 = POP .* y0 .* ( RR_1 - 1 ) ./ RR_1;        
        base_d2 = POP .* y0 .* ( RR_2 - 1 ) ./ RR_2;        
        
        % extrapolate from base_d1 using frac to get base_d for cfn
        base_d = base_d1 + frac * (base_d2 - base_d1);     
        
        % assign base deaths to output array
        base(a,o) = base_d;
        
        % loop through source groups depending on number of variables
        if numVars == 2
            for i = 1:sz(1)
            for j = 1:sz(2)
               % get source contribution
               temp = dJ(i,j);
               % calculate fraction between points for linear extrapolation
               % NOTE: if considering large sources ind_1 and ind_2 should
               % be recalculated w/ (cfn-temp) >= EXP and (cfn-temp) < EXP
               frac  = (cfn - temp - z0)/(zf-z0);    
               % calculate perturbed pert
               pert_d = base_d1 + frac * (base_d2-base_d1);      
               % assign pert d to output array
               D(i,j) = D(i,j) + (base_d-pert_d);                
            end
            end            
        elseif numVars == 3
            for i = 1:sz(1)
            for j = 1:sz(2)
            for k = 1:sz(3)
               % get source contribution
               temp = dJ(i,j,k);
               % calculate fraction between points for linear extrapolation
               % NOTE: if considering large sources ind_1 and ind_2 should
               % be recalculated w/ (cfn-temp) >= EXP and (cfn-temp) < EXP
               frac  = (cfn - temp - z0)/(zf-z0);    
               % calculate perturbed pert
               pert_d = base_d1 + frac * (base_d2-base_d1);      
               % assign pert d to output array
               D(i,j,k) = D(i,j,k) + (base_d-pert_d);                
            end
            end                  
            end
        elseif numVars == 4
            for i = 1:sz(1)
            for j = 1:sz(2)
            for k = 1:sz(3)
            for l = 1:sz(4)
               % get source contribution
               temp = dJ(i,j,k,l);
               % calculate fraction between points for linear extrapolation
               % NOTE: if considering large sources ind_1 and ind_2 should
               % be recalculated w/ (cfn-temp) >= EXP and (cfn-temp) < EXP
               frac  = (cfn - temp - z0)/(zf-z0);    
               % calculate perturbed pert
               pert_d = base_d1 + frac * (base_d2-base_d1);      
               % assign pert d to output array
               D(i,j,k,l) = D(i,j,k,l) + (base_d-pert_d);                
            end
            end                  
            end    
            end
        else
           fprintf(['This function is only set-up to calculate source', ...
                   ' apportionment for 2d, 3d, and 4d data'])
        end 
end
end
end

function [base,D] = o3_health_calc(pop,rate,cfn,dJ,beta,tmrel)

% Function "o3_health_calc" takes inputs of population data (pop),
% mortality rates (rate), a base cost-function value (cfn), and an adjoint 
% contribution array (dJ) and calculates the total and source apportionment of
% heatlh impacts.

% Variables:
%-------------------------------------------------------------------------%
% pop  - m x 1 vector of age-stratified population data
% rate - m x n array of outcome and age seperated mortality rates
% cfn  - scalar cost-function value of o3exposure
% dJ   - 2D, 3D, or, 4D array of sources that contribute to the cost-function 
%        exposure

% Outputs
%-------------------------------------------------------------------------%
% base - total o3 exposure deaths from cost-function
% D    - an array of source apportionment deaths

%% Set-up

% Assignments
amax = length(pop);    % number of age groups
omax = size(rate,2); % number of health outcomes
% get size of source apportionment array
numVars = length(size(dJ)); sz = zeros(numVars,1);
for n = 1:numVars
   sz(n) = size(dJ,n); 
end
% health parameters from GBD 2019
% beta = log(1.06) / 10; tmrel = 32.4;

% zcf     = 32.4; % tmrel from Cohen et al. 2017
% zcf_ub  = 29.1;
% zcf_lb  = 35.7;
% 
% beta    = log(1.06) / 10; % relative risk from Jerret (Cohen et al. supp)
% beta_ub = log(1.10) / 10; % relative risk from Jerret (Cohen et al. supp)
% beta_lb = log(1.03) / 10; % relative risk from Jerret (Cohen et al. supp)

% Initializations
D    = zeros(sz');
base = zeros(amax,omax);

%% Calculate premature deaths

% loop through age groups and outcomes
for a = 1:amax
for o = 1:omax
         
        % get age population
        POP = pop(a);
        % get mortality rate for outcome and age group
        y0 = rate(a,o) / 1E5;
        % get relative risk from log-linear exposure response equation
        if cfn > tmrel; RR = exp(beta * (cfn - tmrel)); else; RR = 1; end
                        
        % calculate cost-function premature deaths for both risks
        base_d = POP .* y0 .* ( RR - 1 ) ./ RR;             
        
        % assign base deaths to output array
        base(a,o) = base_d;
        
        % loop through source groups depending on number of variables
        if numVars == 2
            for i = 1:sz(1)
            for j = 1:sz(2)
               % get source contribution
               temp = dJ(i,j); 
               
               % get relative risk from log-linear exposure response equation
               if cfn > tmrel
                   RR = exp(beta * (cfn - temp - tmrel)); 
               else
                   RR = 1; 
               end         
               
               % calculate cost-function premature deaths for both risks
               pert_d = POP .* y0 .* ( RR - 1 ) ./ RR;             
                       
               % assign pert d to output array
               D(i,j) = D(i,j) + (base_d-pert_d);                
            end
            end            
        elseif numVars == 3
            for i = 1:sz(1)
            for j = 1:sz(2)
            for k = 1:sz(3)
               % get source contribution
               temp = dJ(i,j,k); 
               
               % get relative risk from log-linear exposure response equation
               if cfn > tmrel
                   RR = exp(beta * (cfn - temp - tmrel)); 
               else
                   RR = 1; 
               end         
               
               % calculate cost-function premature deaths for both risks
               pert_d = POP .* y0 .* ( RR - 1 ) ./ RR;  
               
               % assign pert d to output array
               D(i,j,k) = D(i,j,k) + (base_d-pert_d);                
            end
            end                  
            end
        elseif numVars == 4
            for i = 1:sz(1)
            for j = 1:sz(2)
            for k = 1:sz(3)
            for l = 1:sz(4)
               % get source contribution
               temp = dJ(i,j,k,l); 
               
               % get relative risk from log-linear exposure response equation
               if cfn > tmrel
                   RR = exp(beta * (cfn - temp - tmrel)); 
               else
                   RR = 1; 
               end         
               
               % calculate cost-function premature deaths for both risks
               pert_d = POP .* y0 .* ( RR - 1 ) ./ RR;    
               
               % assign pert d to output array
               D(i,j,k,l) = D(i,j,k,l) + (base_d-pert_d);                
            end
            end                  
            end    
            end
        else
           fprintf(['This function is only set-up to calculate source', ...
                   ' apportionment for 2d, 3d, and 4d data'])
        end 
end
end
end

function [base,D] = no_health_calc(pop,rate,cfn,dJ,beta,tmrel)

% Function "o3_health_calc" takes inputs of population data (pop),
% mortality rates (rate), a base cost-function value (cfn), and an adjoint 
% contribution array (dJ) and calculates the total and source apportionment of
% heatlh impacts.

% Variables:
%-------------------------------------------------------------------------%
% pop  - m x 1 vector of age-stratified population data
% rate - m x n array of outcome and age seperated mortality rates
% cfn  - scalar cost-function value of o3exposure
% dJ   - 2D, 3D, or, 4D array of sources that contribute to the cost-function 
%        exposure

% Outputs
%-------------------------------------------------------------------------%
% base - total o3 exposure deaths from cost-function
% D    - an array of source apportionment deaths

%% Set-up

% Assignments
amax = length(pop);    % number of age groups
omax = size(rate,2); % number of health outcomes
% get size of source apportionment array
numVars = length(size(dJ)); sz = zeros(numVars,1);
for n = 1:numVars
   sz(n) = size(dJ,n); 
end

% Initializations
D    = zeros(sz');
base = zeros(amax,omax);

%% Calculate premature deaths

% loop through age groups and outcomes
for a = 1:amax
for o = 1:omax
         
        % get age population
        POP = pop(a);
        % get mortality rate for outcome and age group
        y0 = rate(a,o) / 1E5;
        % get relative risk from log-linear exposure response equation
        if cfn > tmrel; RR = exp(beta * (cfn - tmrel)); else; RR = 1; end
                        
        % calculate cost-function premature deaths for both risks
        base_d = POP .* y0 .* ( RR - 1 ) ./ RR;             
        
        % assign base deaths to output array
        base(a,o) = base_d;
        
        % loop through source groups depending on number of variables
        if numVars == 2
            for i = 1:sz(1)
            for j = 1:sz(2)
               % get source contribution
               temp = dJ(i,j); 
               
               % get relative risk from log-linear exposure response equation
               if cfn > tmrel
                   RR = exp(beta * (cfn - temp - tmrel)); 
               else
                   RR = 1; 
               end         
               
               % calculate cost-function premature deaths for both risks
               pert_d = POP .* y0 .* ( RR - 1 ) ./ RR;             
                       
               % assign pert d to output array
               D(i,j) = D(i,j) + (base_d-pert_d);                
            end
            end            
        elseif numVars == 3
            for i = 1:sz(1)
            for j = 1:sz(2)
            for k = 1:sz(3)
               % get source contribution
               temp = dJ(i,j,k); 
               
               % get relative risk from log-linear exposure response equation
               if cfn > tmrel
                   RR = exp(beta * (cfn - temp - tmrel)); 
               else
                   RR = 1; 
               end         
               
               % calculate cost-function premature deaths for both risks
               pert_d = POP .* y0 .* ( RR - 1 ) ./ RR;  
               
               % assign pert d to output array
               D(i,j,k) = D(i,j,k) + (base_d-pert_d);                
            end
            end                  
            end
        elseif numVars == 4
            for i = 1:sz(1)
            for j = 1:sz(2)
            for k = 1:sz(3)
            for l = 1:sz(4)
               % get source contribution
               temp = dJ(i,j,k,l); 
               
               % get relative risk from log-linear exposure response equation
               if cfn > tmrel
                   RR = exp(beta * (cfn - temp - tmrel)); 
               else
                   RR = 1; 
               end         
               
               % calculate cost-function premature deaths for both risks
               pert_d = POP .* y0 .* ( RR - 1 ) ./ RR;    
               
               % assign pert d to output array
               D(i,j,k,l) = D(i,j,k,l) + (base_d-pert_d);                
            end
            end                  
            end    
            end
        else
           fprintf(['This function is only set-up to calculate source', ...
                   ' apportionment for 2d, 3d, and 4d data'])
        end 
end
end
end




