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

function plotFigure(fP,fA,comodulogram,clabel,dirOut,savename)

    f = figure('visible','off');
    set(f,'DefaultAxesFontName', 'FreeSerif')
    set(f,'DefaultAxesFontSize', 8)
    imagesc(fP,fA,comodulogram)
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) pos(2) pos(3)-0.05 0.98*pos(4)])
    set(gca,'Ydir','normal')
    set(gca,'TickDir','out')
    set(gca,'XTick',fP)
    set(gca,'XTickLabel',fP)
    ylabel('Frequency for amplitude [Hz]')
    xlabel('Frequency for phase [Hz]')
    box off

    pos = get(gca,'Position');
    axes('Position',[0.98*(pos(1)+pos(3)) pos(2) 0 pos(4)]);
    set(gca,'visible','off')
    c = colorbar();
    if sum(sum(comodulogram>0))>0
        caxis([min(min(comodulogram)) max(max(comodulogram))]);
    end
    ylabel(c,clabel)
    cWidth = get(c,'Position');
    cWidth(3) = cWidth(3)*1.1;
    set(c,'Position',cWidth)
    set(c,'TickDir','out')
    box off
        
    print(f, '-dpng', '-r300', [dirOut savename '.png']);
    close(f);
    
    
    
end