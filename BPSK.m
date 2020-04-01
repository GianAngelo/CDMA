function [modData] = BPSK(data)
 % Create a QPSK modulator System object with bits as inputs and Gray-coded signal constellation
%        qpsk_mod = comm.BPSKModulator('BitInput',true);
       Bpsk_mod = comm.BPSKModulator();
 % Modulate
       modData = step(Bpsk_mod,data);
end