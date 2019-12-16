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

eeglab

EEG = pop_loadset('filename', 'TestData.set', 'filepath', './');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
EEG = eeg_checkset( EEG);

chan_fP    = 1;
chan_fA    = 1;
epochs     = 1;
dirOut     = './';
type       = 'continous';
fPstart    = 4;
fPend      = 8;
fPstep     = 1;
fPband     = 1;
fPfiltband = 1;
fAstart    = 20;
fAend      = 120;
fAstep     = 5;
fAfiltband = fPend;
nbCycles   = 3;
nbBins     = 18;
Athresh    = 0.1;
Nboot      = 200;
w          = 5;
peMI       = 95;
pMI        = 95;
pdPAC      = 95;
pPhaseCom  = 95;
artifacts  = true;
plotWithMask    = 1;
plotWithoutMask = 0;

% running Extended Modulation Index Analysis
eMI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPband,fAstart,fAend,fAstep,dirOut,type,nbCycles,nbBins,w,Athresh,Nboot,peMI,...
             pPhaseCom,artifacts,plotWithMask,plotWithoutMask);

% running Modulation Index Analysis
MI_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,nbBins,Nboot,pMI,...
            artifacts,plotWithMask,plotWithoutMask);

% running direct PAC Estmate Analysis
dPAC_calc(EEG,chan_fP,chan_fA,epochs,fPstart,fPend,fPstep,fPfiltband,fAstart,fAend,fAstep,fAfiltband,dirOut,type,Nboot,pMI,...
            artifacts,plotWithMask,plotWithoutMask);
        