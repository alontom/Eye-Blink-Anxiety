% The authors give no warranty for the correct functioning of the software
% and cannot be held legally accountable.

%% Initiate variables

% Enter the number of participants and chunks per participant
num_participants=50; num_chunks=18;

% Choose EEG channels
ch1 = 'Ve1'; ch2 = 'Fp1';

% Set high and low pass filter
hp = 0.5; lp = 20;

% Set minimum peak and minimum distance between consecutive peaks
min_peak = 100; min_distance = 250;

% Enter the number of figures to be plotted
plt=1;

% For each participant collect their blink characteristics
blinks=nan(num_participants,num_chunks); % Total number blinks
blinks_ts=cell(num_participants,num_chunks); % Blink time-series
blink_amp=cell(num_participants,num_chunks); % Blink amplitude
len_in_min=nan(num_participants,num_chunks); % Length in seconds
blink_st=zeros(num_participants,num_chunks,751); % Blink structure
% Blink Rate, Blink Amplitude, Coefficient of Variation (IBI, Amplitude)
[br,ba,brv,bav]=deal(nan(num_participants,num_chunks)); 


%% Detect blinks
for part = 1:num_participants
    % Adding 0 before 1-9
    if part<10
        part_str=['0',num2str(part)];
    else
        part_str=num2str(part);
    end
    % Change to the main data folder
    folder = ['...',part_str];
    % Validate the exsitence of the folder
    if isfolder(folder)
        % Print the participant's number
        disp(part)
        for chunk = 1:num_chunks            
            % Load the EEG dataset (eg. part_01_chunk_1.set)
            file = sprintf('part_%s_chunk_%i.set',part_str,chunk);
            EEG = pop_loadset(file,folder);

            % Compute the length in minutes based on the sampling rate
            len_in_min(part,chunk)=EEG.pnts/(60*EEG.srate);

            % Find the channel indices for 'Ve1' and 'Fp1'
            ch1_index = strcmpi({EEG.chanlocs.labels}, ch1);
            ch2_index = strcmpi({EEG.chanlocs.labels}, ch2);

            % Extract data from 'Ve1' and 'Fp1' channels
            ch1_data = EEG.data(ch1_index, :);
            ch2_data = EEG.data(ch2_index, :);

            % Compute the difference between the two electrodes
            d_data = ch1_data - ch2_data;

            % High and low pass filters         
            d_data_filtered = eegfilt(d_data, EEG.srate, hp,0,0,0,0,'fir1');
            d_data_filtered = eegfilt(d_data_filtered, EEG.srate, 0,lp,0,0,0,'fir1');

            % Find blink time amplitude
            [blink_amp_d,blink_ts_d] = findpeaks(d_data_filtered,"MinPeakHeight",min_peak,"MinPeakDistance",min_distance);
            blink_ts_d = blink_ts_d/EEG.srate; % Convert time to seconds

            % Compute the average blink structure - 250 data points before
            % blink onset up to 500 after
            for b=blink_ts_d
                if b>250 && b<EEG.pnts-500
                    blink_st(part,chunk,:) = squeeze(blink_st(part,chunk,:))'+d_data_filtered(b-250:b+500)/length(blink_ts_d);
                end
            end
           
            % Add up the blinks at the specific chunk
            blinks(part,chunk) = length(blink_ts_d); % Count blinks
            blinks_ts{part,chunk} = blink_ts_d; % Blink time-series
            blink_amp{part,chunk} = blink_amp_d; % Blink amplitude

            % If you want to plot the Fp1-Ve1 filtered with blink regions
            if plt>0
                % Plot data
                figure;
                plot(EEG.times/EEG.srate, d_data_filtered); % XX
                % Paint a 750ms region surrounding the blink peak (250
                % before and 500 after)
                xregion(blink_ts_d - 0.25,blink_ts_d + 0.5, FaceColor='r')
                title(['Participant ',part_str, ch1, ch2]);
                xlabel('Time (s)');
                ylabel('Amplitude (μV)');
                set(gcf,'color','w');
                plt=plt-1;
            end
        end
    end
end
%% Compute measures: BR, BA, BRV, BAV

for part=1:num_participants
    if ~isnan(blinks(part,1))
        for chunk=1:num_chunks
            % BR - blinks per minute
            br(part,chunk)=blinks(part,chunk)/len_in_min(part,chunk);

            % BA - Average blink amplitude
            ba(part,chunk)=mean(blink_amp{part,chunk});

            % Inter-Blink-Interval - Time difference between two consecutive
            % blinks
            ibi=diff(blinks_ts{part,chunk});

            % BRV - SD(IBI)/mean(IBI)
            brv(part,chunk)=std(ibi)/mean(ibi);

            % BAV - SD(Amp.)/mean(Amp.)
            bav(part,chunk)=std(blink_amp_film{part,chunk})/mean(blink_amp_film{part,chunk});
        end
    end
end
