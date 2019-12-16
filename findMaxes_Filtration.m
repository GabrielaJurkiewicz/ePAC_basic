%    Phase-amplitude coupling detection MATLAB plugin
%
%    Copyright (C) 2019 Gabriela Jurkiewicz
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%    Gabriela Jurkiewicz <gabriela.j.jurkiewicz@gmail.com>

function [Maxes,LowFreq] = findMaxes_Filtration(EEG,chan_fP, fP_bins, LimitAmp,order)

    progressbar('on')
    Maxes = {};
    LowFreq = {};
    
    for Idx = 1:size(fP_bins,2)
        
        fP_b = fP_bins(:,Idx)';
        Mx = {};
        LF = {};
        
        for trial = 1:EEG.trials
            
            Fs = EEG.srate;  % Sampling Frequency
            N   = order;  % Order
            Fc1 = fP_b(1);  % First Cutoff Frequency
            Fc2 = fP_b(2);  % Second Cutoff Frequency
            h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
            Hd = design(h, 'butter');
            sigF  = filtfilt(Hd.sosMatrix,Hd.ScaleValues,double(EEG.data(chan_fP,:,trial)));
            
            IdxMax = [];
            if sum((sigF>=LimitAmp)~=0)>0
                [~,IdxMax] = findpeaks(double(sigF),'MINPEAKDISTANCE',floor((0.5/mean(fP_b))*EEG.srate),'MINPEAKHEIGHT', LimitAmp);
            end
            
            Mx{trial} = IdxMax;
            LF{trial} = sigF;
            
            progressbar(((Idx-1)*EEG.trials+trial)/(size(fP_bins,2)*EEG.trials)*100)
            
        end
        
        Maxes{Idx}   = Mx;
        LowFreq{Idx} = LF;
        progressbar(Idx/size(fP_bins,2)*100)
    end
    progressbar('off')
    
end