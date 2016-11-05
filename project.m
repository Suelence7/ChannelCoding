function [] = project()
% Name of Author: George Samir Ibrahim
% GUC ID: 25-10773
% Tutorial: T-05
% Team Members: Only ME
clc;
clear;
close all;
clear all;
% read video file
video = mmreader('.\sv\foreman.avi');
% call stream function applying all prerequsets in the project
% p=0.1 using no channel coding
stream(video, 0.1, 0, 0, 0, 1, '0.1noCC.avi');
% p=0.1 using rate 1/2 convolutional without incremental redundancy
stream(video, 0.1, 1, 0, 0, 1, '0.1noRed.avi');
% p=0.1 using incremental redundancy
stream(video, 0.1, 1, 1, 0, 1, '0.1withRed.avi');
% p=0.001 using no channel coding
stream(video, 0.001, 0, 0, 0, 1, '0.001noCC.avi');
% p=0.001 using rate 1/2 convolutional without incremental redundancy
stream(video, 0.001, 1, 0, 0, 1, '0.001noRed.avi');
% p=0.001 using incremental redundancy
stream(video, 0.001, 1, 1, 0, 1, '0.001withRed.avi');
% Plot of the coded bit error probability against different values of p
% (Assume a range of p between 0.0001 and 0.2) assuming rate 1/2 code 
% given without using incremental redundancy.
stream(video, 0.0001:0.001:0.2, 1, 0, 1, 0, 'withChannelCodingGraph');
% Plot of the coded bit error probability against different values 
% of p (Assume a range of p between 0.0001 and 0.2) using incremental 
% redundancy. And Plot of the throughput against different values 
% of p (Assume a range of p between 0.0001 and 0.2) using incremental
% redundancy.
stream(video, 0.0001:0.001:0.2, 1, 1, 1, 0, 'withIncRedGraph');
end

function cc = stream(video, p, cc, incRed, drawGraphs, saveMovie, output)
% video = video file
% p = probability of error
% incRed = with incremental redundency = 1 oterwise = 0
% cc = with channel coding = 1 otherwise = 0
% drawGraphs = with draw Graphs = 1 otherwise = 0
% saveMovie = to save Movie = 1 otherwise = 0
% output = file name to be saved
% Initialization
trellis = poly2trellis(7,[133 171]); % inialize the trellis with genrator
                                     % code [133 171] with rate 1/2
receivedVideo(1:get(video,'NumberOfFrames')) = struct('cdata', ...
    uint8([]), 'colormap', []); % create movie structure matrix.
% Initializing puncturing sequence
for i = 1:128
    punc89(1,(16*(i-1) + 1):(16*i)) = [1 1 1 0 1 0 1 0 0 1 1 0 1 0 1 0];
    punc45(1,(16*(i-1) + 1):(16*i)) = [1 1 1 0 1 0 1 0 1 1 1 0 1 0 1 0];
    punc23(1,(16*(i-1) + 1):(16*i)) = [1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0];
    punc47(1,(16*(i-1) + 1):(16*i)) = [1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0];
end
n = 1;
% total the avg ratio and throuphput for each video per probabililty
% of error
AvgRatio = zeros(size(p));
avgThroughPut = zeros(size(p));
[fH fW fP] = size(read(video,1));
noOfPackets = (get(video,'NumberOfFrames')*fH*fW*fP*8)/1024;
for pe = p
    % reset the avg ratio and throuphput for each video per probabililty 
    % of error
    ratio = zeros(1,noOfPackets); % avg bitError ratio for each
                                  % RGB for each frame
    throughPut = zeros(1,noOfPackets); % avg throuphput for each RGB for
                                       % each frame 
    noOfPackets = 0;
    n2 = 1;
    % loop on the frames of the video to cut the video into pieces
    for i = 1:get(video,'NumberOfFrames')
        %read each frame
        readFrame = read(video,i);
        % fpages is the RGB colour mixer in 3 integers for the same 
        % pixel in each frame where fH = frame Height, fW = frame width,
        % fP= frame color page.
        [fH fW fP] = size(readFrame);
        % loop on colour RGB in each frame fPages usually here is 3
        for j = 1:fP
            % read color page reshaping the decimal matrix of colour  
            % to row matrix converting decimal representation to binary  
            % in each row reshaping it into on row of bit sequence 
            % prepared to packetizing and being sent
            bCR = reshape(de2bi(reshape((double(readFrame(:,:,j)))',1,...
                (fH*fW)))',1,(fH*fW*8));
            % loop to divide the bit sequence into 1024 bit per packet and
            % sending them to be received by the other side (receiver)
            for k = 1:(ceil((fH*fW*8)/1024))
                if(k*1024 < (fH*fW*8))
                    till = k*1024; % end point in binary color row array
                                   % of the packet to be transmitted
                else
                    till = fH*fW*8;
                end
                % packet is the packet supposed to be transmitted
                packet = bCR((((k-1)*1024) +1):till);
                flag = 1; % flag to make sure if the packet recevied right
                          % or increase error for Incremental Redunency
                error = 0;
                bitCounter = 0; % bit counter is number of bits received 
                                % for each packet
                while flag > 0
                    if (cc < 1)
                        % with no channel coding we just send the packet 
                        % as it is and noise to it
                        rec = bsc(packet,pe);
                        % add number of recevied bits
                        [~,bitCounter] = size(packet);
                        noOfPackets = noOfPackets + 1;
                        flag = 0;
                    else
                        if (incRed < 1)
                            % encode the binary bits using the trellis then 
                            % send frame over channel with channel noise
                            % probabilty equal to zero
                            se = bsc(convenc(packet, trellis),pe);
                            % decode the received binary string using the
                            % trellis
                            rec = vitdec(se,trellis,40,'trunc','hard');
                            % add number of recevied bits
                            [~,bitCounter] = size(se);
                            noOfPackets = noOfPackets + 1;
                            flag = 0;
                        else
                            if(error == 0)
                                se = bsc(convenc(packet, trellis,...
                                    punc89),pe);
                                rec = vitdec(se,trellis,40,'trunc',...
                                    'hard', punc89);
                                % add number of recevied bits
                                [~,bitCounter] = size(se);
                                noOfPackets = noOfPackets + 1;
                            else if(error == 1)
                                    se = bsc(convenc(packet, trellis,...
                                        punc45),pe);
                                    rec = vitdec(se,trellis,40,'trunc',...
                                        'hard', punc45);
                                    % add number of recevied bits
                                    [~,bitCounter] = size(se);
                                    noOfPackets = noOfPackets + 1;
                                else if(error == 2)
                                        se = bsc(convenc(packet,...
                                            trellis,punc23),pe);
                                        rec = vitdec(se,trellis,40,...
                                            'trunc','hard', punc23);
                                        % add number of recevied bits
                                        [~,bitCounter] = size(se);
                                        noOfPackets = noOfPackets + 1;
                                    else if(error == 3)
                                            se = bsc(convenc(packet,...
                                                trellis, punc47),pe);
                                            rec = vitdec(se,trellis,40,...
                                                'trunc','hard', punc47);
                                            % add number of recevied bits
                                            [~,bitCounter] = size(se);
                                            noOfPackets = noOfPackets + 1;
                                        else if(error == 4)
                                                se = bsc(convenc(packet....
                                                    ,trellis),pe);
                                                rec = vitdec(se,trellis,...
                                                    40,'trunc','hard');
                                                % add number of recevied
                                                % bits
                                                [~,bitCounter] = size(se);
                                                noOfPackets = ...
                                                    noOfPackets + 1;
                                            end
                                        end
                                    end
                                end
                            end
                            % calculate the error of the recevied bits
                            e = biterr(rec,bCR((((k-1)*1024) +1):till));
                            % if the packet recevied in high error go for
                            % next punctring
                            if(error < 4)
                                if(e ~= 0)
                                    error = error + 1;
                                    flag = 1;
                                else
                                    flag = 0;
                                end
                            else
                                flag = 0;
                            end
                        end
                    end
                end
                % save the received packet to the frame
                r((((k-1)*1024) +1):till) = rec;
                % average ratio for each packet colors in each frame
                [~,c1]=size(packet);
                ratio(n2) = (biterr(rec,packet)/c1);
                throughPut(n2) = (c1/bitCounter);
                n2 = n2 +1;
            end
            receivedVideo(1,i).cdata(:,:,j) =  reshape(bi2de((reshape(r,...
                8,fH*fW)')),fW,fH)'; % regenerate the frame and assign it
                                     % to it's place in the received video
            fprintf('%s, pe = %f, frame = %f, color = %f\n',output,pe,i,j);
        end
        %average ratio for each frame
    end
    %ratio of error for video per porbabitity of error
    AvgRatio(n)  = (sum(ratio)/noOfPackets);
    avgThroughPut(n) = (sum(throughPut)/noOfPackets);
    n = n + 1;
end
% Save movie
if(saveMovie>0)
movie2avi(receivedVideo,output,'compression', 'None');%Create AVI file
% movieview(output); % play the movie
end
% Draw the graphs and save
if( drawGraphs > 0)
figure(1);
h = plot(p,AvgRatio,'green');
title(output);
xlabel('probability of error');
ylabel('Average bit error ratio');
name = strcat(output,'BitError');
saveas(h,name,'fig');
figure(2);
h = plot(p,avgThroughPut,'green');
title(output);
xlabel('probability of error');
ylabel('Average Throughput');
name = strcat(output,'ThrouphPut');
saveas(h,name,'fig');
end
end