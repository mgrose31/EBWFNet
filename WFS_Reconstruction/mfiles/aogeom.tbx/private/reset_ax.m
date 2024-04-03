function reset_ax
% SYNTAX:
% reset_ax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) MZA Associates Corporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Id: reset_ax.m 3378 2011-03-16 03:38:47Z venet $

%% BEGIN_CODE

dmglobal;

delete(get(gui_main,'currentaxes'));
figure(gui_main);

% We use the Outer Inert Radius (2nd element of 3rd row of matrix) to
% calculate the size of axes

% scaling factors and tic mark increments for x- and y- axis are equal
%%%%%
% (BV, 2011-03-15:) Move axes origin further inward so that tick labels
% don't blank out so frequently (still not ideal, but much improved):
%%% myaxes = axes('position',[0.05,0.05,0.9,0.9],'xlim',[0,1],'ylim',[0,1]);
myaxes = axes('position',[0.1,0.1,0.85,0.85],'xlim',[0,1],'ylim',[0,1]);
%%%%%
axis('equal');
xaxl = get(myaxes, 'xlim');
yaxl = get(myaxes, 'ylim');
xaxlsz = xaxl(2)-xaxl(1);
yaxlsz = yaxl(2)-yaxl(1);

if (op_disp_counts == 0)
        yaxl = [dm_act_radii(3,2)*(-1),dm_act_radii(3,2)];
else
        yaxl = [dm_act_radii(3,2)*(-1),dm_act_radii(3,2) + 0.1];
end

if (xaxlsz >= yaxlsz)
        xaxl = yaxl * (xaxlsz/yaxlsz);
else
        xaxl = [dm_act_radii(3,2)*(-1),dm_act_radii(3,2)];
        yaxl = xaxl * (yaxlsz/xaxlsz);
end
axis('normal');

set(myaxes,'xlim',xaxl,'ylim',yaxl);

if (op_disp_axes)
   set(myaxes,'visible','on');
else
   set(myaxes,'visible','off');
end

set(myaxes, 'DefaultTextFontName', 'times');
set(myaxes, 'DefaultTextFontSize', 12);
set(myaxes, 'DefaultTextFontWeight', 'bold');

% (BV, 2005-07-20:) The following additions seems to cancel a bug (observed using
% Matlab 6.5.1) that caused axes region background to turn white when certain 
% changes (e.g., hide/show changes) were made from the GUI View menu.
% (Change was needed bec. white bgd obscures several graph elements):
set(myaxes, 'Color', [0,0,0]);
set(myaxes, 'XColor', [1,1,1]);
set(myaxes, 'YColor', [1,1,1]);

