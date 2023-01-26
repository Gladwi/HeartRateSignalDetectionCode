basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);
FS = 20000; %samples per second (HZ?)
%timeTotal = 300; %seconds
a = memmapfile([basepath filesep basename '.dat'],'format','int16');
xml = LoadXml([basepath filesep basename '.xml']);


%a = memmapfile([basepath filesep basename '.dat'],[numOfChan, FS*timeTotal],'format','int16');
% Kathryn's version of using the memmmapfile funtion to create an empty
% array/vector that the user can determine the size of the desired
% extracted data



ekgChan = xml.AnatGrps(4).Channels + 1;
%from the xml file in the AnatGrps, it says channel six, becuase neuroscope
%is on a scale that starts from 0, but matlab starts from one, so even
%though it is channel 6 in neuroscope, in matlab that same channel is
%designated as 7. (xml file shows 6, but matlab reads that on [1-infinity
%sacle} so in matlab it is at position 6, while in neuroscope it is
%position 7

numOfChan = 7;
ekg = a.data(ekgChan:numOfChan:1e6); %can use 3e7 or 4.1e6 or 1e6


%ekg = a.data(ekgChan,:);
%Kathryn's version of ekg. If using the Kat_ver of memmapfile use then this
%function is usefule because the data will only contain the channels that
%I am interested in



%the number of channels here signifies that we are pulling only the seventh
%data point for each single timestamp recording. basically 7 variables are
%collected simultaneously and we are only pulling into matlab the 7th one
%the 1e6 is the maximum number of data points we want to pull out (so we
%only want to see up to a million data points

% figure
% plot(ekg);

timestamps = (0:length(ekg)-1)/xml.SampleRate;


%before including the imcluding the timestamp x-axis the data was just
%separated by unit spaces and nothing of value
%The total number of samples we are looking at are 142857 samples
%20000 samples were colleeted every second(20KHz) so we can find how long
%this recording period is by dividing the samples by the sampling rate
%once we do that we now have a scale showing when each sample was collected
%over a 7.14285 second period (use this as the x-axis to collerate each
%data point to its coorect timing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%unfiltered EKG/ECG
%figure 
%plot(timestamps,ekg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%HighPass Filtering The EKG/ECG
doub_ekg = double(ekg);
highfilt_ekg = highpass(doub_ekg, 20, 2e4);
% figure
% plot(timestamps, highfilt_ekg);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LowPass Filtering The EKG/ECG

low_highfilt_ekg = lowpass(highfilt_ekg, 1, 2e4);

% figure
% plot(timestamps, low_highfilt_ekg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%HISTORGRAM
figure
histogram(low_highfilt_ekg, 100);

%this creates a histogram that shows the distrbution of points. this lets
%you see the distribution of points and if there are possible thresholds
%that can be implemented that would better filter your data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Attempting to figure out findpeaks

%[Minima1, MinIdx1] = findpeaks(-low_highfilt_ekg, 20000, 'MinPeakDistance', 0.0700)
%                                                 Hz                   secs
%using minPeakDistance ran into issues with peaks that occured below the
%min distance. It also included minimum peaks that were just noise. I will
%try to use MinPeakHeight instead
[Minima, MinIdx] = findpeaks(-low_highfilt_ekg, xml.SampleRate, 'MinPeakHeight', 2500); 
%%%%%%% figure out how to detect one peak within a certain interval using findpeaks 
%%%%%%% certain interval length should be the theoretical maximum of heart
%%%%%%% rate
Minima = -Minima; %the Minima values come out as positive numbers, so this just sets them as their actual negative numbers

% hold on
figure;
plot(timestamps, low_highfilt_ekg);
% xline(MinIdx); %test to see if index positions line up with minima

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create average R-R interval calculator
%Removes outliers
%Outliers are time differnces less than 0.071 seconds and those more than
%0.194 seconds based on john hopkins univeristy heart rate of 310 - 840 bpm

RR_sum = 0;
true_pos_count = 0;
false_pos_count = 0;
possible_missed_intervals = 0;
RRint_loc_matrix = double.empty(0, 2);
for i = 2:(length(MinIdx) - 1) % the first data point is skipped 
    timeDif =  MinIdx(i+1) - MinIdx(i);
    if timeDif > 0.071 && timeDif < 0.194
        RR_sum = RR_sum + timeDif;
        true_pos_count = true_pos_count+1;
        RRint_loc_matrix(true_pos_count, 1) = (MinIdx(i+1) + MinIdx(i))/2;
        RRint_loc_matrix(true_pos_count, 2) = timeDif;

    elseif timeDif < 0.071
        false_pos_count = false_pos_count + 1;
    else
        possible_missed_intervals = possible_missed_intervals + 1;
    end
end
avg_RR = RR_sum/true_pos_count;

%well now I could have taken the filt_RRs vector and just run it through an
%avg function, but oh well


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BAR GRAPH 
% (time range based on John Hopkins Uni mouse heartrate range 310-840 bpm)
x = categorical({'True Positive', 'Missed Intervals', 'False Positive'});
x = reordercats(x, {'True Positive', 'Missed Intervals', 'False Positive'});
y = [true_pos_count possible_missed_intervals false_pos_count];
bar(x, y)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DataMakeup Percentage
%Outliers are time differnces less than 0.071 seconds and those more than
%0.194 seconds based on john hopkins univeristy heart rate of 310 - 840 bpm
numTimDifTot = true_pos_count + possible_missed_intervals + false_pos_count;
falsePosPerc = (false_pos_count/numTimDifTot) * 100;
possibleMissPerc = (possible_missed_intervals/numTimDifTot) * 100;
truePosPerc = (true_pos_count/numTimDifTot) * 100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create a vector containing the post-processed filtered RR_intervals
%Outliers are time differnces less than 0.071 seconds and those more than
%0.194 seconds based on john hopkins univeristy heart rate of 310 - 840 bpm

j = 1;
filt_RRs = zeros(1,true_pos_count);
for i = 2:(length(MinIdx) - 1)
    timeDif =  MinIdx(i+1) - MinIdx(i);
    if timeDif > 0.071 && timeDif < 0.194
          filt_RRs(j) = timeDif;
          j = j + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SDNN/SDRR while still removing the outliers (in seconds)
%Outliers are time differnces less than 0.071 seconds and those more than
%0.194 seconds based on john hopkins univeristy heart rate of 310 - 840 bpm

SDNN = std(filt_RRs);


%checking that this is working properly

% dif_squarSum = 0;
% for i = 1:length(filt_RRs)
%     dif_squarSum = dif_squarSum + (abs(filt_RRs(i) - avg_RR))^2;
%     SDNNtest = sqrt(dif_squarSum/length(filt_RRs));
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Root Mean Square (in seconds)
rootMeanSq = rms(filt_RRs);

%Test
% squar_filtRRs = filt_RRs.^2;
% sum_squarfiltRRs = 0;
% for i = 1:length(filt_RRs)
%     sum_squarfiltRRs =  sum_squarfiltRRs + squar_filtRRs(i);
% end
% rootMeanSqTest = sqrt(sum_squarfiltRRs/length(filt_RRs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NN50 (code to count how many times the a 50 ms interval exist between
% (cleaned/NN interval) QRS complexes) (exclude first index) (time range based on John Hopkins Uni mouse heartrate
% range)

NN50 = 0;
for i = 2:(length(MinIdx) - 1)
    timeDif = MinIdx(i+1) - MinIdx(i);
    if timeDif > 0.071 && timeDif < 0.194
        if timeDif > 0.050
            NN50 = NN50 + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NN Count Code (time range based on John Hopkins Uni mouse heartrate
% range)
NNcount = 0;
for i = 2:(length(MinIdx) - 1)
    timeDif = MinIdx(i+1) - MinIdx(i);
    if timeDif > 0.071 && timeDif < 0.194
        NNcount = NNcount + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pNN50 Code

pNN50 = NN50 / NNcount;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Method for checking if I am detecting peaks extrememly close to each
%%other

% figure;
% plot(low_highfilt_ekg(1:2e4));
% hold on; plot(MinIdx(1:10)*2e4,Minima(1:10),'ro'); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%addpath(genpath('\\homedir-cifs.nyumc.org\gto219\Personal\MATLAB\buzcode-master'))
%run this code everytime opening matlab so I can access the functions
%created by David
