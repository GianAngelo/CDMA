function [modData] = BPSKdemod(data)
% Create a QPSK modulator System object with bits as outputs and Gray-coded signal constellation
       Bpsk_mod = comm.BPSKDemodulator();
% Demodulate the data using the demodulator object 
       modData = step(Bpsk_mod,data);
end 