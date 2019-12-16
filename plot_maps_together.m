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

function [] = plot_maps_together(x,y,mapa,signal,meanTheta,PC,PhaseBins,savename,cut,plotWithMask)
    
    if plotWithMask
        phase = angle(hilbert(meanTheta));
        mask  = zeros(size(mapa));
        for fa = 1:length(y)
            IDX = zeros(size(phase));
            PH = find(PC(fa,:)>0);
            if ~isempty(PH)
                for ph = 1:length(PH)
                    idxPH = ((phase>PhaseBins(PH(ph))) & (phase<=PhaseBins(PH(ph)+1)));
                    IDX = IDX | idxPH;
                end
            end
            mask(fa,:) = IDX;
        end
    end

    fontNm = 'FreeSerif';
    labelSize = 10;    

    f = figure('visible','off');
    colormap('parula');
    set(f,'DefaultAxesFontName', fontNm)
    set(f,'DefaultAxesFontSize', 10)
    
    subplot(5,1,[1,2,3,4]);
    title('Averaged TF map','FontSize',labelSize,'FontName',fontNm)
    imagesc(x,y,log(mapa))
    if plotWithMask
        if sum(sum(mask))>0
            hold on;
            contour(x,y,mapa.*mask,[eps eps],'.','Color',[0 0 0]/255,'LineWidth',1.5);
        end
    end
    set(gca,'position',[0.1,0.35,0.75,0.60])
    set(gca,'Ydir','normal')
    set(gca,'XTick',[])
    set(gca,'TickDir','out')
    box off
    ylabel('Frequency [Hz]','FontSize',labelSize,'FontName',fontNm)

    c = colorbar();
    set(c,'position',[0.87,0.35,0.03,0.60])
    set(c,'TickDir','out')
    if sum(sum(mapa))>0
        caxis([min(min(log(mapa))) max(max(log(mapa)))])
    end
    box off
    ylabel(c,'log(power)','FontSize',labelSize,'FontName',fontNm)

    subplot(5,1,5);
    title('Averaged high- and low-freq signal','FontSize',labelSize,'FontName',fontNm) 
    plot(x,signal,'Color',[1,184,202]/255,'LineWidth',2)
    hold on
    plot(x,meanTheta,'Color',[0,0,0]/255,'LineWidth',1)
    set(gca,'position',[0.1,0.11,0.75,0.17])
    set(gca,'XTick',round(linspace(-cut,cut,7),2))
    set(gca,'TickDir','out')
    box off
    xlim([-cut, cut])
    xlabel('Time [s]','FontSize',labelSize,'FontName',fontNm)
    ylabel('Amplitude','FontSize',labelSize,'FontName',fontNm)
    
    print(f, '-dpng', '-r300', [savename '.png']);
    close(f);
    
end