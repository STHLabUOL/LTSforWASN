% Coh = defineCoherenceArbitr_func(type,P,Dist,c,fs,nfft)
%   generate matrix with desired spatial coherence
%   implemented for use with Habets' coherent noise generator
%   (arbitrary geometry, does not need to be an array)
%
% Input: 
%   type            'spherical' or 'cylindrical'
%   P               # of microphones
%   Dist            distance matrix for microphones [m]
%   c               speed of sound [m/s]
%   fs              sampling frequency [Hz]
%   nfft             FFT-length
%
% Output:
%   Coh             P x P x nfft/2+1 spatial coherence matrix
%
% written by Philipp Thuene
% modified August 23, 2018

function Coh = defineCoherenceArbitr_func(type,P,Dist,c,fs,nfft)

    w = 2*pi*fs*(0:nfft/2)/nfft;
    Coh = zeros(P,P,nfft/2+1);
    for p = 1:P
        for q = 1:P
            if p == q
                Coh(p,q,:) = ones(1,1,nfft/2+1);
            else
                if strcmpi(type,'spherical')
                    Coh(p,q,:) = sinc(w*Dist(p,q)/(c*pi));

                elseif strcmpi(type,'cylindrical')
                    Coh(p,q,:) = besselj(0,w*Dist(p,q)/c);

                else
                    error('Unknown noise field.')
                end
            end
        end
    end
    
end