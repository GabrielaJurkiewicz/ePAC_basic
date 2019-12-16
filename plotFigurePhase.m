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

function plotFigurePhase(fP,fA,PHASE,PhaseComod,clabel,dirOut,savename)

    PhaseComod = reshape(PhaseComod,[length(fA),length(fP)*length(PHASE)]);

    % rysowanie i podstawowe ustawienia osi
    f = figure('visible','off');
    set(f,'DefaultAxesFontName', 'FreeSerif')
    set(f,'DefaultAxesFontSize', 8)
    imagesc(1:1:length(fP)*length(PHASE),fA,PhaseComod)
    xlim([0.5 (length(fP)*length(PHASE))+1])
    set(gca, 'Ydir', 'normal')
    set(gca,'TickDir','out')
    
    % zwezanie mapki zeby colorbar sie zmiescil
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) pos(2) pos(3)-0.05 0.98*pos(4)])
    
    % rysowanie pionowych bialych linii rozdzielacjacych fP
    for l = 2:length(fP)
        hold on 
        line([(l-1)*length(PHASE)+0.5 (l-1)*length(PHASE)+0.5],get(gca,'YLim'),'LineWidth',1.5,'Color',[134 118 220]/255)
    end
    
    % ustawienia osi X i Y mapki; X faza wyswietlana na gorze mapy
    phase_ticks = repmat(PHASE',length(fP),1)';
    phase_labels = {};
    idxP = find((phase_ticks == -pi) | (phase_ticks == 0));
    for pt = 1:length(idxP)
        if phase_ticks(idxP(pt))==-pi
            if pt==1
                phase_labels = [phase_labels '-\pi'];
            else
                phase_labels = [phase_labels '\pi'];
            end
        else
            phase_labels = [phase_labels '0'];
        end
    end
    idxP = [idxP length(fP)*length(PHASE)+1];
    phase_labels = [phase_labels '\pi'];
    set(gca,'XTick',idxP-0.5)
    set(gca,'XTickLabel',phase_labels)
    set(gca,'xaxisLocation','top')
    box off
    xlabel('Phase')
    ylabel('Frequency for amplitude [Hz]')
    
    % ustawienia colorbar
    pos = get(gca,'Position');
    axes('Position',[0.98*(pos(1)+pos(3)) pos(2) 0 pos(4)]);
    set(gca,'visible','off')
    c = colorbar();
    if sum(sum(PhaseComod))>0
        caxis([min(min(PhaseComod)) max(max(PhaseComod))]);
    end
    ylabel(c,clabel)
    cWidth = get(c,'Position');
    cWidth(3) = cWidth(3)*1.1;
    set(c,'Position',cWidth)
    set(c,'TickDir','out')
    box off
    
    % dodatkowa os na dole mapy zawierajaca podpisy czestosci fP
    fPlabel = {};
    for p = 1:length(fP)
        fPlabel = [fPlabel num2str(fP(p))];
    end
    axes('Position',[pos(1) pos(2) pos(3) 0]);
    set(gca,'XTick',((length(PHASE)/(length(PHASE)*length(fP)):length(PHASE)/(length(PHASE)*length(fP)):1)-length(PHASE)/(length(PHASE)*length(fP)*2)))
    set(gca,'XTickLabel',fPlabel)
    set(gca,'TickDir','out')
    set(gca,'xaxisLocation','bottom')
    set(gca,'TickLength',[0 0])
    xlabel('Frequency for phase [Hz]')
    box off
    
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    set(findall(gcf,'-property','FontName'),'FontName','FreeSerif')
    
    print(f, '-dpng', '-r300', [dirOut savename '.png']);
    close(f);
    
end