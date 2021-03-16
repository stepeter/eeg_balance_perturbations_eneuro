%%EEG processing
eeg_files = dir('data_release_final/*.set');
ica_struct_files = dir('data_release_ICA_structs/*ICA_struct.mat');

for i=1:length(eeg_files)
    % Load data and ICA_STRUCT
    EEG = pop_loadset('filename', eeg_files(i).name, 'filepath','data_release_final/');
    load(['data_release_ICA_structs/' ica_struct_files(i).name]);
    
    % Remove bad channels using ICA_STRUCT
    EEG=pop_select(EEG,'channel',ICA_STRUCT.good_chans); %remove bad channels
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
    pop_saveset(EEG, 'filename', eeg_files(i).name, 'filepath', 'data_release_denoised/','savemode','twofiles');
end
