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
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%    Gabriela Jurkiewicz <gabriela.j.jurkiewicz@gmail.com>

function [Maps, MeanSignal, MessMaps, MeanFPsignal] = calcMeanMap(signal, mx, fPsignal, Fs, cut, w, fAstart, fAend, fAstep, fA, fP_b, Nboot, plotWithMask)

    %% ---------------- PREPARE STRUCTURES ---------------------
    num = 0;
    Maps         = zeros(size(fA,2),floor(cut*Fs)*2+1);
    MessMaps     = zeros(Nboot,size(fA,2),floor(cut*Fs)*2+1);
    MeanSignal   = zeros(1,floor(cut*Fs)*2+1);
    MeanFPsignal = zeros(1,floor(cut*Fs)*2+1);

    %% ---------------- CALC MORLET WAVELET TF MAP ---------------------
    [map,~,~] = tf_cwt(hilbert(signal),fAstart,fAend+fAstep,Fs,w,fAstep,0);

    %% ------------- ADDING MAPS ALIGNED TO MAXIMA -------------- 
    while ~isempty(mx)
        
        idx  = mx(1,1);
        num = num + 1;
        mx(mx<(idx+floor(2*cut*Fs))) = [];
        
        Maps = Maps + map(:,idx-floor(cut*Fs):idx+floor(cut*Fs));
        MeanSignal = MeanSignal + signal(idx-floor(cut*Fs):idx+floor(cut*Fs));
        MeanFPsignal  = MeanFPsignal + fPsignal(idx-floor(cut*Fs):idx+floor(cut*Fs));

        if plotWithMask
            for nb = 1:Nboot
                id = floor((rand()-1.0/2)*Fs/mean(fP_b));
                scale = (rand(1)-0.5)*0.1/0.5; % form -0.1 to 0.1
                scale = (1+scale);
                idxBoot = floor(idx*scale) + id;
                MAX = floor(size(map,2)*scale);
                while (idxBoot-floor(cut*Fs)<=0 || idxBoot+floor(cut*Fs) > MAX)
                    id = floor((rand()-1.0/2)*Fs/mean(fP_b));
                    scale = (rand(1)-0.5)*0.1/0.5; % form -0.1 to 0.1
                    scale = (1+scale);
                    idxBoot = floor(idx*scale) + id;
                    MAX = floor(size(map,2)*scale);
                end
                tmp = imresize(map,[size(map,1) floor(size(map,2)*scale)]);
                MessMaps(nb,:,:) = squeeze(MessMaps(nb,:,:)) + tmp(:,idxBoot-floor(cut*Fs):idxBoot+floor(cut*Fs));
            end
        end
    end

    %% ------------- MINIMAL NUMBER OF TRIALS CONDITION -------------------
    if (num < 3)
        num = 0;
        Maps         = zeros(size(fA,2),floor(cut*Fs)*2+1);
        MessMaps     = zeros(Nboot,size(fA,2),floor(cut*Fs)*2+1);
        MeanSignal   = zeros(1,floor(cut*Fs)*2+1);
        MeanFPsignal = zeros(1,floor(cut*Fs)*2+1);
    end

    %% ------------- NORMALIZATION -------------------
    if (num > 0)
        MeanSignal   = MeanSignal/num;
        MeanFPsignal = MeanFPsignal/num;
        Maps         = Maps/num;
        if plotWithMask
            MessMaps     = MessMaps/num;
        end
    end
    
end
