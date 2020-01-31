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

function [] = run_continuous_eMI(EEG,Epochs,HighFreqSignal,Maxes,LowFreq,fP_bins,fP,fAstart,fAend,fAstep,margines,dirOut,nbCycles,w,nbBins,Nboot,pPhaseCom,peMI,ID,plotWithMask,plotWithoutMask)

    fA = fAstart:fAstep:fAend;
    PhaseBins = zeros(1,nbBins);
    dPhi      = 2*pi/nbBins;
    for b = 1:nbBins 
        PhaseBins(b) = -pi+(b-1)*dPhi; 
    end
    
    for epoch = 1:length(Epochs)
 
        disp(['Epoch ' num2str(Epochs(epoch)) ' / ' num2str(EEG.trials)])
        progressbar('on')
        Comod = zeros(length(fA),size(fP_bins,2));
        ComodSurr = zeros(Nboot,length(fA),size(fP_bins,2));
        PhaseComod = zeros(length(fA),nbBins,size(fP_bins,2));

        results = matfile([dirOut 'v' ID '_eMI_epoch' num2str(Epochs(epoch)) '_results.mat'],'Writable',true);
        results.fA = fA;
        results.fP = fP;
        results.fP_bins = fP_bins;
        
        for i = 1:size(fP_bins,2)

            %% ---------------- ADJUST DATA TO SPECIFIC fP  ---------------------
            fP_b  = fP_bins(:,i)';
            cut   = nbCycles/(2*fP(1,i));
            time  = linspace(-cut,cut,floor(cut*EEG.srate)*2+1);

            signal   = squeeze(HighFreqSignal(1,:,Epochs(epoch)));
            fPsignal = LowFreq{i}{1,Epochs(epoch)};
            mxs      = Maxes{i}{1,Epochs(epoch)};
            mx       = mxs((mxs>(floor((cut+margines)*EEG.srate)))&(mxs<(floor(size(signal,2))-floor((cut+margines)*EEG.srate))));

            %% ---------------- CALCULATE TF MAPS + MI -------------------
            [Map, MeanSignal, MessMap, MeanLowFreq] = calcMeanMap(signal, mx, fPsignal, EEG.srate, cut, w, fAstart, fAend, fAstep, fA, fP_b, Nboot, plotWithMask);
            [MI,MISurr,PC] = calc_ModulationIndex(Map,MessMap,MeanLowFreq,nbBins,PhaseBins,Nboot,pPhaseCom,plotWithMask);
            clear MessMap
            Comod(:,i) = MI;
            ComodSurr(:,:,i) = MISurr;
            PhaseComod(:,:,i) = PC;          
            clear MI MISurr PC
            
            %% -------------------- SAVE AND PLOT RESULTS + CLEAN -----------------------------
            eval(['results' num2str(fP(1,i)) 'Hz.time = time;'])
            eval(['results' num2str(fP(1,i)) 'Hz.Map = Map;'])
            eval(['results' num2str(fP(1,i)) 'Hz.MeanSignal = MeanSignal;'])
            eval(['results' num2str(fP(1,i)) 'Hz.MeanLowFreq = MeanLowFreq;'])
            eval(['results.results' num2str(fP(1,i)) 'Hz = results' num2str(fP(1,i)) 'Hz;'])
            eval(['clear results' num2str(fP(1,i)) 'Hz'])
            clear Map MeanSignal MeanLowFreq
            
            progressbar(i/size(fP_bins,2)*100)
        end
        
        progressbar('off')
        disp('----------------------- eMI: ploting and saving results -----------------------------------------------------')
        
        %% ---------------- eMI STATISTIC --------------------------------------
        if plotWithMask
            MX  = max(max(ComodSurr,[],2),[],3);
            clear ComodSurr
            tresh = prctile(MX,peMI);
            mask = (Comod>=tresh);

            comodulogramStat = zeros(size(Comod));
            comodulogramStat(mask) = Comod(mask);

            pval = zeros(size(Comod));
            PhaseComod_afterStat = zeros(size(PhaseComod));
            for fp = 1:length(fP)
                for fa = 1:length(fA)
                    pval(fa,fp) = sum(MX>=Comod(fa,fp))/Nboot;
                    if mask(fa,fp)==1
                        PhaseComod_afterStat(fa,:,fp) = PhaseComod(fa,:,fp);
                    end
                end
            end
        end
        
        %% ---------------- SAVE AND PLOT RESULTS + CLEAN --------------------------------------
        if plotWithoutMask
            plotFigure(fP,fA,Comod,'eMI',dirOut,['v' ID '_eMI_epoch' num2str(Epochs(epoch)) '_comodulogram'])
        end
        results.comodulogram = Comod;
        clear Comod
        
        if plotWithMask
            plotFigure(fP,fA,comodulogramStat,'eMI',dirOut,['v' ID '_eMI_epoch' num2str(Epochs(epoch)) '_comodulogramStat'])
            plotFigurePhase(fP,fA,PhaseBins,PhaseComod,'eMI',dirOut,['v' ID '_eMI_epoch' num2str(Epochs(epoch)) '_phaseComodulogram']) 
            plotFigurePhase(fP,fA,PhaseBins,PhaseComod_afterStat,'eMI',dirOut,['v' ID '_eMI_epoch' num2str(Epochs(epoch)) '_phaseComodulogramStat']) 
            results.phaseComodulogram = PhaseComod;
            clear PhaseComod
            results.phaseComodulogramStat = PhaseComod_afterStat;
            clear PhaseComod_afterStat
            results.comodulogramStat = comodulogramStat;
            clear ComodulogramMaxStat
            results.pval = pval;
            clear pval
        end

        results.phase_vector = PhaseBins; 
        matfile([dirOut 'v' ID '_eMI_epoch' num2str(Epochs(epoch)) '_results.mat'],'Writable',false);
        clear results

        load([dirOut 'v' ID '_eMI_epoch' num2str(Epochs(epoch)) '_results.mat']);
        for i = 1:size(fP_bins,2)
            cut   = nbCycles/(2*fP(1,i));
            time  = linspace(-cut,cut,floor(cut*EEG.srate)*2+1);
            eval(['time = results' num2str(fP(1,i)) 'Hz.time;'])
            eval(['Map = results' num2str(fP(1,i)) 'Hz.Map;'])
            eval(['MeanSignal = results' num2str(fP(1,i)) 'Hz.MeanSignal;'])
            eval(['MeanLowFreq = results' num2str(fP(1,i)) 'Hz.MeanLowFreq;'])
            if plotWithMask
                eval(['PhaseComod_afterStat = results' num2str(fP(1,i)) 'Hz.MeanLowFreq;'])
            else
                phaseComodulogramStat = zeros(length(fA),nbBins,size(fP_bins,2));
            end
            plot_maps_together(time,fA,Map,MeanSignal,MeanLowFreq,squeeze(phaseComodulogramStat(:,:,i)),[phase_vector pi],...
                                 [dirOut 'v' ID '_eMI_epoch' num2str(Epochs(epoch)) '_results' num2str(fP(1,i)) 'Hz'],cut,plotWithMask)
        end
        
    end
   
    if length(Epochs)>1
        plotContinuousSummary(EEG,fP,fA,dirOut,ID,(100-peMI)/100,plotWithMask,plotWithoutMask)
    end
    
end
