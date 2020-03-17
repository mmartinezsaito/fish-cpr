% spm_help QUESTION
% see DEM and FieldMap toolboxes
% spm_jobman('interactive',matlabbatch)
% visualization: spm_image, spm_orthviews, spm_spm;;
% contrast validity: spm_conman, spm_SpUtil, spm_sp(942)
% orthogonalization: spm_orth (Gram-Schmidt orthogona(norma)lization);


%clearvars -except <> 
close all %#ok<CLSCR>;

load ~/Yandex.Disk/CPR/sazan_fitted_cfl.mat

spmpath = '~/spm12'; % ~/spm12 or Toolboxes/spm8 or ~/Data/Applications/spm12
addpath(spmpath); 
%spm fmri
spm('defaults', 'FMRI');
global default_params; 
default_params = spm_get_defaults;
spm_jobman('initcfg');

%% setting parameters

% T2*-weighted imaging protocol
TR = 2.28;
nslices = 40; sliceMatrix = [76 76];
FoV = 228; % mm

% utils
comp = @(s) s(:,1) - mean(s(:,2:3), 2);
sust = @(s) -abs(6 -s(:,1) -s(:,2) -s(:,3)); 

% base folders
bdir = '~/Data/CPR/NIfTI-1/'; 
anadir = '~/Yandex.Disk/CPR/FmriAnalysis';      
jobdir = [anadir '/jobs']; 
behlogsdir = '~/Yandex.Disk/CPR/logs';
nifti1dir = '~/Data/CPR/NIfTI-1/'; % anadir
%nifti1dir = '~/Data/CPR/vasily_NIfTI-1/'; 

%addpath(genpath(bdir))

% is social model first orthogoonalized contrast?
mt = {'social', 'nonsoc'};
issoc1 = 0;  % 1, 0
isrpe = 0;
if issoc1, sstrtail = 'sn';
  modtyp1 = mt{1}, modtyp2 = mt{2}; rtyp1 = comp; rtyp2 = sust;
else,      sstrtail = 'ns';
  modtyp1 = mt{2}, modtyp2 = mt{1}; rtyp1 = sust; rtyp2 = comp;
end
if ~isrpe, sstrtail = [sstrtail 'r'], else sstrtail, end

% Current analysis directory
tilde2dir = @(x) regexprep(x, '~', '/home/mario'); %
%   classical,bayesian; canonical,informed,fir; categorical,parametric
cad = ['Classical_CanonicalHrf_Pmod']; 
antyp = [tilde2dir(anadir) filesep cad filesep];
polex = 1; % polynomial order of modulating parameter (order of polynomial expansion where 0 is none)

usefm = 0; % use fieldmap?

runjob.import = 0;
getnifti1info = 0;

% Preprocessing
runjob.slicetimingcorrection         = 0;
runjob.realignment.estimation        = 0;
runjob.realignment.estimation_unwarp = 0;
runjob.coregistration                = 0;
runjob.segmentation                  = 0;
runjob.normalization                 = 0;
runjob.smoothing                     = 0;

% First-level model specification and estimation
runjob.specify_1st.classical         = 0;
runjob.estimate_1st.classical        = 0;
runjob.estimate_1st.bayesian         = 0; 
runjob.contrasts                     = 0;

% Second-level model specification and estimation
runjob.specify_2nd.classical         = 1;
runjob.estimate_2nd.classical        = 1;
runjob.estimate_2nd.bayesian         = 0; 
 
% Bayesian Model Selection?
runbms = 0;
% Compute LogEvMatrix
if runbms, Compute_LogEvMatrix_BMS, end

% jobs better done manually
runjob.results = 0;

% subject ids
sids = 2:51;  
soc_sids = [2:10 20 21 26 27 30 31 34 35 37 38 44 45 47 49 51];
nsoc_sids = setdiff(sids, soc_sids);
discarded_sid = [30 43 46]; % discarded 
multcol_sid = [30]; % discarded because hamper analysis with multicollinearity: s30 has identical parameters, variables, and regressors for s,n models
corrupted_data = [43 46]; 
culled_sids = setdiff(sids, discarded_sid)


%% import
if runjob.import
    load([jobdir '/import_50']);
    output_list = spm_jobman('run', matlabbatch);
end


%% Loop over subjects while running prespecified jobs
starttime = clock, startday = date

for i =  culled_sids  
    i
    dsstr = sprintf('sub%02u', i);
    asstr = sprintf(['sub%02u_' sstrtail], i);
    
    fdir = tilde2dir([nifti1dir dsstr '/fMRI']);
    sdir = tilde2dir(fullfile(nifti1dir, dsstr, '/sMRI'));
    fmdir = tilde2dir(fullfile(nifti1dir, dsstr, '/fMapping'));
    fimg = dir([fdir '/f*.img']); fhdr = dir([fdir '/f*.hdr']);
    simg = dir([sdir '/s*.img']); shdr = dir([sdir '/s*.hdr']);
    fmimg = dir([fmdir '/s*.img']); fmhdr = dir([fmdir '/s*.hdr']);
        
    % Get header info
    if getnifti1info
        Vf = spm_vol(tilde2dir([fdir filesep fimg(1).name]));
        Vs = spm_vol(tilde2dir([sdir filesep simg(1).name]));
        Vfm = spm_vol(tilde2dir([fmdir filesep fmimg(1).name]));
        Vf.dim, Vs.dim, Vfm.dim
    end
            
    % Slice timing correction
    if runjob.slicetimingcorrection
        clear matlabbatch SPM files
        %filepath = spm_select('FPList', fdir, '^f.*img');
        for j=1:length(fimg)
            files{j} = [fdir filesep fimg(j).name];
        end
        load([jobdir '/slicetiming_i']);
        matlabbatch{1}.spm.temporal.st.scans{1} = files';
        matlabbatch{1}.spm.temporal.st.nslices = nslices; 
        matlabbatch{1}.spm.temporal.st.tr = TR; 
        matlabbatch{1}.spm.temporal.st.ta = TR - TR/nslices; 
        matlabbatch{1}.spm.temporal.st.so = 1:40; % ascending
        matlabbatch{1}.spm.temporal.st.refslice = round(nslices/2); 
        output_list = spm_jobman('run', matlabbatch);
    end    
    
    % Realignment
    if runjob.realignment.estimation
        clear matlabbatch SPM files;
        prefix = 'a';
        for j=1:length(fimg)
            files{j} = [fdir filesep prefix fimg(j).name];
        end        
        load([jobdir '/realign_i']);
        %if ispc, files = files; end
        matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = files';
        output_list = spm_jobman('run', matlabbatch);
        copyfile(output_list{1}.sess.rpfile{1}, [sdir '/movepars.txt']);
        rupf = matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix;
    end
    if runjob.realignment.estimation_unwarp
        clear matlabbatch SPM files
        % Create field map and unwarp EPIs
        s = [homedir '/Data/s09_eokorokova/Okorokova_epi/DICOM/PA000000/ST000000/SE000000/IM000000']; % alternative to NIfTI-1
        Dfm = dicominfo(s);
        % type FieldMap.man; Fieldmap_ngui.m
        % Defaults defined in pm_defaults.m
        % Loading field map data
        global pm_def;
        pm_def.INPUT_DATA_FORMAT = 'PM'; % the input must be either 1 or 2
          % phase and magnitude image pairs. The units for the phase
          % images MUST BE RADIANS BETWEEN +pi and -pi       hm  
        RE = regexp(Vgp{i,1}.descrip, 'T?=(.*?)ms\/', 'tokens');
        refmt1 = regexp(Vgp{i,1}.descrip, 'TE=(.*?)ms\/', 'tokens');
        refmt2 = regexp(Vgp{i,2}.descrip, 'TE=(.*?)ms\/', 'tokens');
        te1 = str2double(refmt1{1});  % 5.19ms
        te2 = str2double(refmt2{1});  % 7.65ms
        pm_def.SHORT_ECHO_TIME = te1;
        pm_def.LONG_ECHO_TIME = te2;
        % Unwrapping
        pm_def.UNWRAPPING_METHOD = 'Mark3D';
        pm_def.PAD = 0;
        pm_def.WS = 1;
        % Saving the field map: fpm_NAME-OF-FIRST-INPUT-IMAGE.img
        % Or loading a precalculated field map: fpm_*img
        % Convert the field map to a VDM (Voxel Displacement Map)
        retrot = regexp(Vf{i}.descrip, 'TE=(.*?)ms\/', 'tokens');
        tert = 1 / Dfm.Private_0019_1028 * 100;  % retrot * Vf{i}.dim(1);    % 32.6ms, "Bandwidth Per Pixel Phase Encode"
        pm_def.TOTAL_EPI_READOUT_TIME = tert;
        pm_def.K_SPACE_TRAVERSAL_BLIP_DIR = -1;
        pm_def.EPI_BASED_FIELDMAPS = 1;
        % Create field map 
        %{
        % Field map preprocess: function for using the FM in batch scripts
        addpath ~/spm12/toolbox/FieldMap
        epi_dir = fdir;
        epifm = 1;
        kdir = -1;
        mask = 0;
        match = 1;
        VDM = FieldMap_preprocess({gmdir, gpdir}, epi_dir, ...
            [te1, te2, epifm, tert, kdir, mask, match]);
        %}        
        % Phase and magnitude data
        load([jobdir '/fm_mp_i']);
        mbfm.defaults.defaultsval.et = [te1 te2];
        mbfm.defaults.defaultsval.maskbrain = 0; 
        mbfm.defaults.defaultsval.blipdir = -1;
        mbfm.defaults.defaultsval.tert = tert;
        mbfm.defaults.defaultsval.epifm = 0;
        mbfm.defaults.defaultsval.ajm = 0;
        mbfm.defaults.defaultsval.uflags.method = 'Mark3D';
        mbfm.defaults.defaultsval.uflags.fwhm = 10;
        mbfm.defaults.defaultsval.uflags.pad = 0;
        mbfm.defaults.defaultsval.uflags.ws = 1;
        mbfm.defaults.defaultsval.mflags.template = cellstr([spmpath filesep 'FieldMap/T1.nii']);
        mbfm.defaults.defaultsval.mflags.fwhm = 5;
        mbfm.defaults.defaultsval.mflags.nerode = 2;
        mbfm.defaults.defaultsval.mflags.ndilate = 4;
        mbfm.defaults.defaultsval.mflags.thresh = 0.5;
        mbfm.defaults.defaultsval.mflags.reg = 0.02;
        mbfm.matchvdm = 1;
        mbfm.sessname = 'session';
        mbfm.writeunwarped = 1;
        mbfm.matchanat = 0;    
        mbfm.anat = ''; 
        mbfm.shortphase = {spm_select('FPList', gpdir, '^s2015.*1\.img')};
        mbfm.shortmag = {spm_select('FPList', gmdir, '^s2015.*\.img')};
        mbfm.longphase = {spm_select('FPList', gpdir, '^s2015.*2\.img')};
        mbfm.longmag = {spm_select('FPList', gmdir, '^s2015.*\.img')};
        mbfm.session.epi = {[fdir filesep fimg(1).name]};
        
        matlabbatch{1}.spm.tools.fieldmap.phasemag.subj = mbfm;
        fm_output_list = spm_jobman('run', matlabbatch);     
                
        % Realign, and unwarp using VDM
        prefix = 'a';
        for j=1:length(fimg)
            files{j} = [fdir filesep prefix fimg(j).name];
        end        
        load([jobdir '/realign_unwarp_i']);
        %if ispc, files = files'; end
        matlabbatch{1}.spm.spatial.realignunwarp.data.scans = files';
        matlabbatch{1}.spm.spatial.realignunwarp.data.pmscan = {spm_select('FPList', gpdir, '^vdm5_scs[0-9].*\.img')};
        output_list = spm_jobman('run', matlabbatch);        
        copyfile(output_list{1}.sess.rpfile{1}, [sdir '/movepars.txt']);
        rupf = matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix;
    end
    if runjob.realignment.estimation || runjob.realignment.estimation_unwarp % plot realignment parameters
        rp = spm_load([sdir '/movepars.txt']); %select the rp*.txt file
        figure
        subplot(2,1,1), plot(rp(:,1:3)), set(gca,'xlim',[0 size(rp,1)+1])
        title([dsstr ' translation']), ylabel('mm'), legend(gca, {'right' 'frontal' 'dorsal'})
        subplot(2,1,2), plot(rp(:,4:6)), set(gca,'xlim',[0 size(rp,1)+1])
        title([dsstr ' rotation']), ylabel('°'), legend(gca, {'pitch' 'roll' 'yaw'})
        print([tilde2dir(anadir) filesep 'rp_' dsstr], '-dpsc')       
    end
    
    % Coregistration
    if runjob.coregistration
        clear matlabbatch SPM files
        load([jobdir filesep 'coreg_i']);
        matlabbatch{1}.spm.spatial.coreg.estimate.ref = {spm_select('FPList', fdir, '^mean.*\.img$')};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {[sdir filesep simg.name]};
        output_list = spm_jobman('run', matlabbatch);
    end
    
    % Segmentation
    if runjob.segmentation
        clear matlabbatch SPM files
        load([jobdir filesep 'segment_i']);
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {[sdir filesep simg.name]};
        output_list = spm_jobman('run',matlabbatch);
    end
    
    % Normalization
    if runjob.normalization
        clear matlabbatch SPM files
        prefix = 'ra';
        if usefm, prefix(1) = 'u'; end
        for j=1:length(fimg)
            files{j} = [fdir filesep prefix fimg(j).name];            
        end
        load([jobdir '/normalize_i']);
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {spm_select('FPList', sdir, '^y_.*')};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = files';
        matlabbatch{2}.spm.spatial.normalise.write.subj.def = {spm_select('FPList', sdir, '^y_.*')};
        matlabbatch{2}.spm.spatial.normalise.write.subj.resample = {spm_select('FPList', sdir, '^ms.*nii$')};
        output_list = spm_jobman('run', matlabbatch);
    end
    
    % Smoothing
    if runjob.smoothing
        clear matlabbatch SPM files
        prefix = 'wra';
        if usefm, prefix(2) = 'u'; end
        for j=1:length(fimg)
            files{j} = tilde2dir([fdir filesep prefix fimg(j).name]);
        end
        load([jobdir '/smooth_i']);
        matlabbatch{1}.spm.spatial.smooth.data = files';
        output_list = spm_jobman('run', matlabbatch);
    end
        
    
    % Specification: 1st level
    if runjob.specify_1st.classical
        clear matlabbatch SPM files
        prefix = 'swra'; 
        if usefm, prefix(3) = 'u'; end            
        for j=1:length(fimg)
            files{1, j} = [fdir filesep prefix fimg(j).name];
        end        
        load([jobdir '/specify1_i_sazan']);    
        if exist([antyp asstr]) ~= 7, mkdir(antyp, asstr); end        
        matlabbatch{1}.spm.stats.fmri_spec.dir = {[antyp asstr]};
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = files';
        
        dmnrows = size(T{T.sid==i,:}, 1);
        dmnd_col = sub2ind([dmnrows 3], [1:dmnrows]', T{T.sid==i,'demanded'}); 
        
        mbcond = matlabbatch{1}.spm.stats.fmri_spec.sess.cond;
        mbcond(1).name = 'Choice';  
        mbcond(1).onset = T{T.sid==i, 'dem_time'};
        mbcond(1).duration = 0;                    
        mbcond(1).pmod(1) = struct('name', {['AnticipatedValue_' modtyp1]}, ...
          'param', {zscore(T{T.sid==i, ['val_' modtyp1]}(dmnd_col))}, 'poly', {polex});
        mbcond(1).pmod(2) = struct('name', {['AnticipatedValue_' modtyp2]}, ...
          'param', {zscore(T{T.sid==i, ['val_' modtyp2]}(dmnd_col))}, 'poly', {polex});
        mbcond(2).name = 'Reaping';  
        mbcond(2).onset = T{T.sid==i, 'rea_time'};
        mbcond(2).duration = 0;            
        mbcond(2).pmod(1) = struct('name', {'Reward'}, ...
          'param', {zscore(T{T.sid==i, 'reaped'})}, 'poly', {polex});
        mbcond(3).name = 'Shrinkage';  
        mbcond(3).onset = T{T.sid==i, 'shr_time'};
        mbcond(3).duration = 0;
        if isrpe
          mbcond(3).pmod(1) = struct('name', {['RPE_' modtyp1]}, ...
            'param', {zscore(T{T.sid==i, ['rpe_' modtyp1]}(dmnd_col))}, 'poly', {polex});
          mbcond(3).pmod(2) = struct('name', {['RPE_' modtyp2]}, ...
            'param', {zscore(T{T.sid==i, ['rpe_' modtyp2]}(dmnd_col))}, 'poly', {polex});
        else
          mbcond(3).pmod(1) = struct('name', {['R_' modtyp1]}, ...
            'param', {zscore(rtyp1(T{T.sid==i, {'reaped' 'shrink1' 'shrink2'}}))}, 'poly', {polex});
          mbcond(3).pmod(2) = struct('name', {['R_' modtyp2]}, ...
            'param', {zscore(rtyp2(T{T.sid==i, {'reaped' 'shrink1' 'shrink2'}}))}, 'poly', {polex});
        end
        mbcond(3).pmod(3) = struct('name', {'CprValueLoss'}, ...
          'param', {zscore(sum(T{T.sid==i, {'shrink1' 'shrink2'}}, 2))}, 'poly', {polex});
        mbcond(4).name = 'CprRegrowthOrDepletion';  
        mbcond(4).onset = T{T.sid==i, 'end_time'};
        mbcond(4).duration = 0;            
        mbcond(4).pmod(1) = struct('name', {'StockSize'}, ...
          'param', {zscore(T{T.sid==i, 'stock'})}, 'poly', {polex});           
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = mbcond; 
        
        % Head motion covariates
        moveparstruc = dir([sdir '/movepars.txt']); %([fdir '/rp*txt']);
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[moveparstruc.folder filesep moveparstruc.name]};
        
        % HRF derivatives
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];  
        %keyboard; spm_jobman('interactive', matlabbatch);
        output_list = spm_jobman('run', matlabbatch);
    end
    
    % Estimation: 1st level
    if runjob.estimate_1st.classical
        clear matlabbatch SPM;
        file = tilde2dir([antyp asstr '/SPM.mat']);
        load([jobdir '/estimate_i']);
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {file};
        output_list = spm_jobman('run', matlabbatch);
    elseif runjob.estimate_1st.bayesian
        clear matlabbatch SPM;
        file = tilde2dir([antyp filesep asstr '/SPM.mat']);
        
        load([jobdir '/estimate_i']); 
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 0;
        %{
        %matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.space.volume.block_type = 'Slices';
        mask = {tilde2dir([anadir0 '/RoiMasks/Striatum_PPC.nii,1'])}; 
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.space.clusters.mask = mask;
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.space.clusters.block_type = 'Subvolumes';
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.signal = 'UGL'; % unweighted graph-laplacian 
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.ARP = 3;
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.noise.UGL = 1;
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.LogEv = 'Yes';        % LogEvidence
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.anova.first = 'No';
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.anova.second = 'Yes';
        matlabbatch{1}.spm.stats.fmri_est.method.Bayesian.gcon = struct('name', {}, 'convec', {});       
        %}        
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {file};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                
        %save(tilde2dir([anadir filesep folder filesep bdir(i).name '/VBbatch' num2str(i) '.mat']),'matlabbatch')
        output_list = spm_jobman('run', matlabbatch);
    end
    
    % Contrasts
    if runjob.contrasts
        clear matlabbatch SPM;
        file = tilde2dir([antyp asstr '/SPM.mat']);
        load([jobdir '/contrast_i']);
        matlabbatch{1}.spm.stats.con.spmmat = {file};

        mbconsess = matlabbatch{1}.spm.stats.con.consess;  
        mbconsess{1}.tcon.name =  'Choice';
        mbconsess{1}.tcon.convec = [1 0 0 0 0 0 0 0 0 0 0];
        mbconsess{2}.tcon.name = ['ChoicePmod: AnticipatedValue ' modtyp1];
        mbconsess{2}.tcon.convec = [0 1 0 0 0 0 0 0 0 0 0];
        mbconsess{3}.tcon.name = ['ChoicePmod: AnticipatedValue ' modtyp2];
        if i==30, mbconsess{3}.tcon.convec = [0 1 0 0 0 0 0 0 0 0 0];          
        else,     mbconsess{3}.tcon.convec = [0 0 1 0 0 0 0 0 0 0 0];
        end
        mbconsess{4}.tcon.name =  'Reaping';
        mbconsess{4}.tcon.convec = [0 0 0 1 0 0 0 0 0 0 0];
        mbconsess{5}.tcon.name =  'ReapingPmod: Reward'; 
        mbconsess{5}.tcon.convec = [0 0 0 0 1 0 0 0 0 0 0];
        mbconsess{6}.tcon.name =  'Shrinkage'; 
        mbconsess{6}.tcon.convec = [0 0 0 0 0 1 0 0 0 0 0];
        mbconsess{7}.tcon.name =  ['ShrinkagePmod1: ' modtyp1];
        mbconsess{7}.tcon.convec = [0 0 0 0 0 0 1 0 0 0 0];
        mbconsess{8}.tcon.name =  ['ShrinkagePmod2: ' modtyp2];
        mbconsess{8}.tcon.convec = [0 0 0 0 0 0 0 1 0 0 0];
        mbconsess{9}.tcon.name =  'ShrinkagePmod: CprValueLoss';
        if i==30, mbconsess{9}.tcon.convec = [0 0 0 0 0 0 0 1 0 0 0];
        else,     mbconsess{9}.tcon.convec = [0 0 0 0 0 0 0 0 1 0 0];
        end
        mbconsess{10}.tcon.name = 'CprRegrowthOrDepletion';
        mbconsess{10}.tcon.convec =[0 0 0 0 0 0 0 0 0 1 0];
        mbconsess{11}.tcon.name = 'CprRegrowthOrDepletionPmod: StockSize';
        mbconsess{11}.tcon.convec =[0 0 0 0 0 0 0 0 0 0 1];
        matlabbatch{1}.spm.stats.con.consess = mbconsess;
                
        % F-tests
        load(file);
        %  SPM.xCon is added ONLY after running contrasts!
        %  SPM.xCon contains  0 SPM-generated contvec/mats: 
        %                    11 defined-by-me contvec/mats
        %                    11 contvec/mats in total
        %  contvec{i} = SPM.xCon{i}.c;
        % if a contrast matrix is outside the row space of the design matrix, run fails
        
        output_list = spm_jobman('run', matlabbatch);
        
        %{
        % Diagnosing contrasts validity
        size(SPM.xX.X)
        figure, imagesc(corr(SPM.xX.X)); colorbar
        figure, imagesc(abs(corr(SPM.xX.X))); colorbar
        surfc(corr(SPM.xX.X)); colorbar, shading interp % surfl
        surfc(corr(abs(SPM.xX.X))); colorbar, shading interp 
        meshc(corr(SPM.xX.X)); colorbar, shading interp
        imagesc(abs(corr(SPM.xX.X))>0.9 & corr(SPM.xX.X)<1); colorbar
        [i j] = find(corr(SPM.xX.X)<1 & corr(SPM.xX.X)>0.99); [i j]'
        imagesc(corr(SPM.xX.xKXs.X)); colorbar       
        assert(rank([SPM.xX.xKXs.X; cvLotW_L]) == rank(SPM.xX.xKXs.X));
        assert(rank([SPM.xX.X; cvLotW_L]) == rank(SPM.xX.X));
        % parametric modulators can render contrasts invalid
        %}
        figure, imagesc(abs(corr(SPM.xX.X))), colorbar, title(dsstr)
        print(gcf, '-depsc', [antyp asstr filesep 'cm'])
                
% When specifying each condition as a separate parametric modulator in a 
% rapid event-related design, the resulting regressors become highly 
% correlated (in the extreme case, they become linearly-dependent). 
% This means that the automatic serial orthogonalization causes later 
% modulators to become closer to zero. This can even produce a column in 
% the design matrix that is zero, and hence the model becomes inestimable 
% (and an error results).
% The solution in such cases is to turn off the (hard-coded) serial 
% orthogonalization of parametric modulators. This involves commenting out 
% the following lines in two separate spm functions: 
%  line  229     in spm_get_ons.m, and 
%  lines 285-287 in spm_fMRI_design.m
% If your resulting modulators are linearly-dependent, you will not be able 
% to estimate certain contrasts (namely those that don't sum to zero) - but 
% this doesn't matter if you are always interested in *differences* between
% conditions, rather than the unique effect of each.        
    end
end


% Specification: 2nd level
if runjob.specify_2nd.classical
    clear matlabbatch SPM;  
    culled = ~ismember(sids, discarded_sid);
    soc_culled = ismember(sids, setdiff(soc_sids, discarded_sid));
    nsoc_culled = ismember(sids, setdiff(nsoc_sids, discarded_sid));
            
    conimg = dir([antyp 'sub02_sn/con_*']);
                
    % Create dirs if non-existent
    for j=1:length(conimg) % osTtest        
        condir_soc = ['c' sscanf(conimg(j).name, 'con_%4s.nii') '_' sstrtail '_soc_' num2str(sum(soc_culled))];
        condir_nsoc = ['c' sscanf(conimg(j).name, 'con_%4s.nii') '_' sstrtail '_nsoc_' num2str(sum(nsoc_culled))];
        if exist([antyp condir_soc], 'file') ~= 7, mkdir(antyp, condir_soc); end
        if exist([antyp condir_nsoc], 'file') ~= 7, mkdir(antyp, condir_nsoc); end        
    end
    for j=1:length(conimg)    % owa    
        condir = ['c' sscanf(conimg(j).name, 'con_%4s.nii') '_' sstrtail '_owa'];
        if exist([antyp condir], 'file') ~= 7, mkdir(antyp, condir); end
    end
    if 0 %cjaimn = {'con_0050.nii', 'con_0051.nii'};
        condir = 'j0001';
        if exist([fpath condir], 'file') ~= 7, mkdir(fpath, condir); end
    end
            
    % OneSampleTtest
    for j=1:length(conimg)
        condir{1} = ['c' sscanf(conimg(j).name, 'con_%4s.nii') '_' sstrtail '_soc_' num2str(sum(soc_culled))];
        condir{2} = ['c' sscanf(conimg(j).name, 'con_%4s.nii') '_' sstrtail '_nsoc_' num2str(sum(nsoc_culled))];
                
        i = 0; files = [];
        for s = setdiff(soc_sids, discarded_sid), i = i + 1;          
            files{1}{i} = [antyp sprintf('sub%02u',s) '_' sstrtail filesep conimg(j).name];
        end
        i = 0;
        for s = setdiff(nsoc_sids, discarded_sid), i = i + 1;                    
            files{2}{i} = [antyp sprintf('sub%02u',s) '_' sstrtail filesep conimg(j).name];
        end
        
        for i = 1:2
          load([jobdir '/specify2_1sTtest_i']);
          matlabbatch{1}.spm.stats.factorial_design.dir = {[antyp condir{i}]};
          matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = files{i}';
          output_list = spm_jobman('run', matlabbatch);
        end
    end
    
    % OneWayAnova
    for j=1:length(conimg)
      condir = ['c' sscanf(conimg(j).name, 'con_%4s.nii') '_' sstrtail '_owa'];
      i1 = 0; i2 = 0; files = [];
      for s = setdiff(sids, discarded_sid)
        if     ismember(s, soc_sids), i1 = i1 + 1;          
          files{1}{i1} = [antyp sprintf('sub%02u',s) '_' sstrtail filesep conimg(j).name];
        elseif ismember(s, nsoc_sids), i2 = i2 + 1;          
          files{2}{i2} = [antyp sprintf('sub%02u',s) '_' sstrtail filesep conimg(j).name];  
        end
      end
        
      load([jobdir '/specify2_1wAnova_i']);
      matlabbatch{1}.spm.stats.factorial_design.dir = {[antyp condir]};
      matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = files{1}';
      matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = files{2}';
      matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0; % measurements assumed independent between levels
      matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1; % measurements in each level assumed to have unequal variance
      output_list = spm_jobman('run', matlabbatch);
    end
       
    % Conjunction analysis (through OneWayAnova)
    %{
    load([jobdir '/specify2_owAnova_i']);
    condir = 'j0001';
    i = 1;
    for s = culled_loop
        files1{i} = [fpath filesep sprintf('sub%02u',s) filesep cjaimn{1}];
        files2{i} = [fpath filesep sprintf('sub%02u',s) filesep cjaimn{2}];
        i = i + 1;
    end
    matlabbatch{1}.spm.stats.factorial_design.dir = {[fpath condir]};
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = files1';
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = files2';
    matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0; % 0: independent
    matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1; % 1: unequal variance
    output_list = spm_jobman('run', matlabbatch);    
    %}
end

% Estimation: 2nd level
if runjob.estimate_2nd.classical || runjob.estimate_2nd.bayesian
    clear matlabbatch SPM;
    a = dir(antyp);    
    i = 1; while i <= length(a)
        ro = regexp(a(i).name, '^([a-z][0-9]{4})_([a-z]*)', 'tokens');
        if ~a(i).isdir || isempty(ro) || ~strcmp(ro{1}{2}, sstrtail)
            a(i) = []; 
        else, i = i + 1;
        end        
    end
    for i=1:length(a)            
        load([jobdir '/estimate_i']);
        file = tilde2dir([antyp a(i).name '/SPM.mat']);
        if runjob.estimate_2nd.classical
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;        
        elseif runjob.estimate_2nd.bayesian
            matlabbatch{1}.spm.stats.fmri_est.method.Bayesian = 1;
        end
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {file};
        output_list = spm_jobman('run', matlabbatch);
    end
end


% Print time summary
endtime = clock, endday = date
fprintf('\n\n\t RUNNING UPTIME')
fprintf('\n\t ============== ')
fprintf(strcat(['\n\tStart: ', startday, ' at ', num2str(starttime(4)), ':', num2str(starttime(5))]))
fprintf(strcat(['\n\tEnd:   ', endday, ' at ', num2str(endtime(4)), ':', num2str(endtime(5))]))
fprintf('\n')

%spm_image('init','niftifile.nii') % display a nifti-file
%spm_orthviews('Reposition', [22 0 -21]) % jump to a certain coordinate
%fs = get(0,'Children');
%res = getframe(fs(1));
%imwrite(res.cdata, 'niftifile.jpg'); % save the window as JPG
% spm_regions extract timeseries
% spm_graph extract and plot timeseries

