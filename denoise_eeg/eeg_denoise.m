%% EEG denoising pipelien
% Start eeglab
eeglab_pth = '.../eeglab13_5_4b/'; % EEGLAB directory
root_pth = 'BIDS/'; % top-level data directory
if ~exist('ALLCOM')
    addpath('plot_fxns/');
    PLTFUNCS.start_eeglab(eeglab_pth)
end

n_eeg_chans = 128; % number of EEG electrodes (same for all participants)
eeg_files = dir([root_pth '*/*/*/sub*_ses-01_task*.set']);
n_subjs = length(eeg_files); clear eeg_files;

for i=1:n_subjs
    eeg_files = dir([root_pth 'sub-' num2str(i,'%03.f') '/*/*/sub-' ...
                     num2str(i,'%03.f') '*_ses*_task*.set']);
    for j=1:length(eeg_files)
        % Load data and ICA_STRUCT
        EEG = pop_loadset('filename', eeg_files(i).name,...
                          'filepath',eeg_files(j).folder);

        % Remove bad channels
        EEG=pop_select(EEG,'channel',...
                       EEG.etc.good_chans(EEG.etc.good_chans<=n_eeg_chans));
        EEG_orig = EEG;

        %Denoising Methods:
        %1) ASR to remove ungrounding events (use st. dev. of 20)
        EEGtemp=clean_asr(EEG,20); EEG.data=EEGtemp.data; EEGtemp=[];

        %2) Run EEMD and CCA to remove large components in 1st IMF
        [EEG,IMF1]=cca_EEMDrem_1stIMF(EEG);

        %Remove externals (last 5 channels)
        EEG=pop_select(EEG,'nochannel',(EEG.nbchan-4):EEG.nbchan);

        %Re-reference and interpolate
        EEG=pop_reref(EEG,[]);
        EEG=pop_interp(EEG,EEG_orig.chanlocs,'spherical');

        % Save dataset
        pop_saveset(EEG, 'filename', eeg_files(i).name, 'filepath',...
                    'BIDS_denoised/','savemode','twofiles');
    end
end
