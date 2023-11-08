function RIR = createRIR_funcAC(room,spkrs,mics,T60,c,len,fs)

    % function RIR = createRIR_func(room,spkrs,mics,T60,c,len,fs)
    %
    % Creates room impulse responses (RIR) using a modified image source
    % method.
    %
    % Requires Lehmann's 'new' toolbox in MATLAB's path variable.
    %
    % Input:
    %   room    room dimensions [3 x 1] in m
    %   spkrs   speaker positions [3 x Q] for Q speakers in m
    %   mics    microphone positions [3 x P] for P microphones in m
    %   T60     reverberation time in s
    %   c       speed of sound in m/s
    %   len     length of RIRs in taps
    %   fs      sampling frequency in Hz
    %
    % Output:
    %   RIR     {Q,P} cell array containing RIRs [len x 1] for all speaker
    %           / microphone combinations
    %
    % References:
    %   Allen, J. B. & Berkley, D. A. 
    %   Image method for efficiently simulating small room acoustics 
    %   Journal of the Acoustical Society of America, 1979, vol. 65, no. 4
    %
    %   Lehmann, E. A. & Johansson, A. M. 
    %   Prediction of energy decay in room impulse responses simulated with an image-source model 
    %   Journal of the Acoustical Society of America, 2008, vol. 124, no. 1
    %
    % written by Philipp Thuene
    % modified March 9, 2012
    
    Q = length(spkrs(1,:));
    P = length(mics(1,:));
    RIR = cell(Q,P);

    % absorption / reflection coefficients
    alpha = gensig.ISM_AbsCoeff('T60',T60,room,ones(1,6),'LehmannJohansson','c',c);
%     alpha = ISM_AbsCoeff('T60',T60,room,ones(1,6),'MillingtonSette','c',c);
%                           % the reflections are stronger that the direct acoustic path
%     beta = (1-alpha).^1.5;  % for T60>0.6s
%     beta = (1-alpha);  % for T60>0.4s
    beta = sqrt(1-alpha);   % already for T60>0.2s
    
    for p = 1:P 
        for q =1:Q    
            % get RIR
            tic
%             display(['From source ',int2str(q),'/',int2str(Q),' to sink ',int2str(p),'/',int2str(P),':']);
            RIR{q,p} = gensig.ISM_RoomResp(fs,beta,'T60',T60,spkrs(:,q),mics(:,p),room,'c',c,'MaxDelay',len/fs);
%             display(['Time spent on RIR calculation: ' num2str(toc), 's.'])
%             display(' ')
        end
    end

end