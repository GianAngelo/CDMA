%CDMA Decoding
%Gian Angelo Tria
%ECE 408: Wireless Communications
clear;
clear all;
clc;
%% Filter and Donwnsample
load('Rcvd_Tria .mat');

%Root Raise Cosine Filter
%Given
B_RCOS = [.0038;.0052;-.0044;-.0121;-.0023;.0143;.0044;-.0385;-.0563;... 
            .0363;.2554;.4968;.6025;.4968;.2554;.0363;-.0563;-.0385;... 
            .0044;.0143;-.0023;-.0121;-.0044;.0052;.0038];
filtered = filter(B_RCOS,1,Rcvd); %%Root Raised Cosine filter used to 
%Downsampling
downsampled = downsample(filtered,4);
%% M Sequence Generation
taps = [1 1 1 0 0 0 0 1];
m_flip_shift = m_seq(taps);

%% Finding Starting Index
indexing=filter(fliplr(reshape([1-2*m_flip_shift;zeros(3,length(m_flip_shift))],1,[])),1,Rcvd);
figure;
mag_indexing = abs(indexing);
plot(mag_indexing)
title('Correlation Between Msequence and Signal')

%% Applying PN Sequencce
downsampled_new = Rcvd(1032:4:end);
repeat_m = repmat(1-2*m_flip_shift,1,ceil(length(downsampled_new)/255));
post_pn = downsampled_new.*repeat_m(1:length(downsampled_new)); 
scatterplot(post_pn);
title('Signal after applying PN Sequence')
%% Fixing Frequency and Phase Shift
angle_post_pn = angle(post_pn);
rotator = cos(-angle_post_pn) + 1j*sin(-angle_post_pn);
fixed = rotator.*post_pn; % At this point everything is on real axis
scatterplot(fixed);
title('Frequency Fixed')
%% Walsh Channel Orthogonal Spreading
h=hadamard(8);
% Find number of complete frames and extract it
copies=(floor(length(fixed)./255));
fixed_new=fixed(1:(copies*255));
Reshape_1=reshape(fixed_new,255,[]);
Reshape_2=reshape(Reshape_1(1:192,:),[],8);
% Each column has 8 chips
Reshape_3=reshape(Reshape_2,8,[]);
decoded = Reshape_3.'*h;
figure
imagesc(abs(decoded))
title('Decoded Signal After Walsh')
demod = pskdemod(decoded,2);
%% Decoding
%Binary to Decimal and then Decimal to Character
binary = transpose(demod(:,6));
for i = 1:length(binary)/8
    characters(i) = bi2de(binary((i-1)*8+1:i*8),'right-msb');
end
decoded_message = char(characters)

%% Functions
function [m_sequence] = m_seq(poly)
    shift = zeros(1, length(poly));
    shift(end) = 1;
    m_sequence = zeros(1,255);
    m_sequence(end) = shift(length(shift));
    for i = 1:length(m_sequence)-1
        shift_end = shift(length(shift));
        for j = length(shift):-1:2
            shift(j) = mod(shift(j-1) + poly(j)*shift_end, 2);
        end
        shift(1) = shift_end;
        m_sequence(length(m_sequence) - i) = shift(end);
    end
end
