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

function [] = plotContinousSummary(EEG,fP,fA,dirOut,ID,alpha,methodName,plotWithMask,plotWithoutMask)

    directoryOut = [dirOut '/v' ID '_' methodName '_SUMMARY/'];
    [status,message,~] = mkdir(directoryOut);
    
    for fp = 1:length(fP)
        num = 0;
        ll = dir(dirOut);
        ll = {ll.name};
        
        for epoch  = 1:EEG.trials
            if strcmp(methodName,'eMI')
                if sum(strcmp(ll,[dirOut 'v' ID '_' methodName '_epoch' num2str(epoch) '_results.mat']))>0
                    load([dirOut 'v' ID '_' methodName '_epoch' num2str(epoch) '_results.mat'])
                    if num == 0
                        MIxtime = zeros(EEG.trials,length(fA));
                        MIpvalxtime = zeros(EEG.trials,length(fA));
                        MIStatxtime = zeros(EEG.trials,length(fA));
                    end
                    fpidx = find(fP==fP(fp));
                    if plotWithoutMask
                        MIxtime(epoch,:) = comodulogram(:,fpidx)';
                    end
                    if plotWithMask
                        MIpvalxtime(epoch,:) = pval(:,fpidx)';
                        MIStatxtime(epoch,:) = comodulogramStat(:,fpidx)';
                    end
                    num = 1;
                end
            else
                if sum(strcmp(ll,['Epoch' num2str(epoch)]))>0
                    load([dirOut 'Epoch' num2str(epoch) '/' 'comodulogram.mat'])
                    load([dirOut 'Epoch' num2str(epoch) '/' 'comodulogramStat.mat'])
                    load([dirOut 'Epoch' num2str(epoch) '/' 'pval.mat'])
                    load([dirOut 'Epoch' num2str(epoch) '/' 'fP.mat'])
                    load([dirOut 'Epoch' num2str(epoch) '/' 'fA.mat'])
                    if num == 0
                        MIxtime = zeros(EEG.trials,length(fA));
                        MIpvalxtime = zeros(EEG.trials,length(fA));
                        MIStatxtime = zeros(EEG.trials,length(fA));
                    end
                    fpidx = find(fP==fP(fp));
                    if plotWithoutMask
                        MIxtime(epoch,:) = comodulogram(:,fpidx)';
                    end
                    if plotWithMask
                        MIpvalxtime(epoch,:) = pval(:,fpidx)';
                        MIStatxtime(epoch,:) = comodulogramStat(:,fpidx)';
                    end
                    num = 1;
                end
            end
        end

        Time = EEG.pnts/EEG.srate*(1:1:EEG.trials);
        units = 's';
        if plotWithMask
            [~, FDRmask] = fdr(MIpvalxtime,alpha,'nonParametric');
            afterFDR = (MIxtime.*FDRmask)';
        end
        if plotWithoutMask
            noFDR = MIxtime';
        end

        %% ------------------- PLOTING AND SAVING ----------------------
        if (EEG.pnts/EEG.srate*EEG.trials>60)
            Time = Time/60;
            units = 'min';
        end
        
        Titles = {};
        if plotWithMask
            Titles = [Titles 'Stat'];
        end
        if plotWithoutMask
            Titles = [Titles ''];
        end
        for i = 1:length(Titles)
            title = Titles{i};
            switch title
                case ''
                    img = noFDR;
                case 'Stat'
                    img = afterFDR;
            end
            % rysowanie i podstawowe ustawienia osi
            f = figure('visible','off');
            set(f,'DefaultAxesFontName', 'FreeSerif')
            set(f,'DefaultAxesFontSize', 8)
            imagesc(Time,fA,img)
            set(gca, 'Ydir', 'normal')
            set(gca,'TickDir','out')
            box off
            xlabel(['Time [' units ']'])
            ylabel('Frequency for amplitude [Hz]')

            % zwezanie mapki zeby colorbar sie zmiescil
            pos = get(gca,'Position');
            set(gca,'Position',[pos(1) pos(2) pos(3)-0.05 0.98*pos(4)])

            % ustawienia colorbar
            pos = get(gca,'Position');
            axes('Position',[0.98*(pos(1)+pos(3)) pos(2) 0 pos(4)]);
            set(gca,'visible','off')
            c = colorbar();
            if sum(sum(img>0))>0
                caxis([min(min(img)) max(max(img))]);
            end
            cWidth = get(c,'Position');
            cWidth(3) = cWidth(3)*1.1;
            set(c,'Position',cWidth)
            set(c,'TickDir','out')
            box off

            set(findall(gcf,'-property','FontSize'),'FontSize',8)
            set(findall(gcf,'-property','FontName'),'FontName','FreeSerif')

            print(f, '-dpng', '-r300', [directoryOut num2str(fP(fp)) 'Hz_' title '.png']);
            close(f);
        end
    end

end