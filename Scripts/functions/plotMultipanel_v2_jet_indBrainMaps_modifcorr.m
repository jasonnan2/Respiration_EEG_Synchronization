function fig = plotMultipanel_v2_jet_indBrainMaps_modifcorr(hm, X, colorlim,labelnames)

[hm.fvLeft,hm.fvRight, hm.leftH, hm.rightH] = geometricTools.splitBrainHemispheres(hm.cortex);

%X = 10000*X; %modify the data range

fig = figure('Color',[1 1 1]);
set(0,'DefaultFigureWindowStyle' , 'normal')
fig.Position(3:4) = [642   755];
% % mx = prctile((X(:)), 95);
% % mn = prctile((X(:)), 5);
% % if mx <= 0
%     mx = max(X(:)); %prctile((X(:)), 97.5);
%     mn = min(X(:)); %prctile((X(:)), 2.5);
% % end
% if mx == mn
if colorlim
    mx = 100; %prctile(abs(X(:)), 80); % 1; %50;
    mn = -100; %-prctile(abs(X(:)), 80); %1; %-50;
   
else
    
    mx = prctile(abs(X(:)), 60)*1.5;
    mn = -mx; %prctile(abs(X(:)), 10);
%     mx = max(X(:)); %prctile((X(:)), 99.5);
%     mn = min(X(:)); %prctile((X(:)), .5);
%     mx = max([abs(mn) mx]);
%     mn = -mx;
    if mn == mx && mn == 0
        mx = max(abs(X(:)));
        mn = -mx;
        if mn == mx && mn == 0
            mn = -0.05;
            mx = 0.05;
        end
    end
end
% end
% L = 10;
% indexValue = 0;     % value for which to set a particular color
% topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
% indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
% bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
% % Calculate where proportionally indexValue lies between minimum and
% % maximum values
% largest = (max(X(:)));
% smallest = (min(X(:)));
% index = L*(indexValue-smallest)/(largest-smallest);
% % Create color map ranging from bottom color to index color
% % Multipling number of points by 100 adds more resolution
% customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
%             linspace(bottomcolor(2),indexColor(2),100*index)',...
%             linspace(bottomcolor(3),indexColor(3),100*index)'];
%         % Create color map ranging from index color to top color
% % Multipling number of points by 100 adds more resolution
% customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
%             linspace(indexColor(2),topColor(2),100*(L-index))',...
%             linspace(indexColor(3),topColor(3),100*(L-index))'];
% customCMap = [customCMap1;customCMap2];  % Combine colormaps

AxesH = axes('Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off', ...
    'YLimMode', 'manual', 'YLim',  [0, 1], ...
    'XTick',    [],       'YTick', [], ...
    'NextPlot', 'add', ...
    'HitTest',  'off');


const = 2;
ncols = size(X,2)+const;
%%
for i = 1:ncols-const

ax = axes('Position', [ (1/ncols)*i*1-0.1    0.8200    0.3000    0.12]);% subplot(511);
patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X(:,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([-90 90])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))
if i == ncols-const
%cb=colorbar;
%cb.Position = cb.Position + 1e-10;
end
% 
%xlabel('Superior','FontWeight','bold')
% ax.Position([3 4]) = [0.82 0.25];
% ax.Position([1 3 4]) = [0.11 0.40 0.3037];
% ax.Position = [0.1500    0.6300    0.3000    0.2000];

% title('top view');

ax = axes('Position', [ (1/ncols)*i*1-0.06    0.75    0.22    0.05]);
patch('vertices',hm.fvLeft.vertices,'faces',hm.fvLeft.faces,'FaceVertexCData',X(hm.leftH,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([-180 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))

%title('Left')
%zlabel('Lateral','FontWeight','bold')
% ax.Position([1 3 4]) = [0.1 0.4091 0.0837];
% ax.Position = [0.1500    0.5500    0.3000    0.0500];
% ax.Position([3 4]) = [0.75 0.06];

% 
ax = axes('Position',[ (1/ncols)*i*1-0.06    0.65    0.22    0.05]);
patch('vertices',hm.fvRight.vertices,'faces',hm.fvRight.faces,'FaceVertexCData',X(hm.rightH,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);

set(ax,'Clim',[mn mx]);
view([0 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))

% ax.Position([1 3 4]) = [0.1 0.3991 0.0837];
% ax.Position = [0.1500    0.4000    0.3000    0.0500];
% ax.Position([3 4]) = [0.75 0.06];

%title('Right')

ax = axes('Position',[ (1/ncols)*i*1-0.07    0.55   0.22    0.055]);
patch('vertices',hm.fvLeft.vertices,'faces',hm.fvLeft.faces,'FaceVertexCData',X(hm.leftH,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([0 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))
% 
% zlabel('Medial','FontWeight','bold')

ax = axes('Position',[ (1/ncols)*i*1-0.07    0.45    0.22    0.055]);
patch('vertices',hm.fvRight.vertices,'faces',hm.fvRight.faces,'FaceVertexCData',X(hm.rightH,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([-180 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');
colormap(bipolar(512, 0.99))


set(findall(fig,'type','axes'),'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','off', 'visible','off')
% fig.Position(3:4) = [405   755];


text((1/ncols)*i*1,0.97,labelnames{i},'units','normalized','Parent', AxesH);


end


text(0.05,0.85,'top','units','normalized','Parent', AxesH);
text(0.05,0.73,'left','units','normalized','Parent', AxesH);
text(0.05,0.57,'right','units','normalized','Parent', AxesH);

text(0.02,0.2,'medial-left','units','normalized','Parent', AxesH);
text(0.02,0.07,'medial-right','units','normalized','Parent', AxesH);


% figure
% 
% ax = axes('Position',[ (1/ncols)*i*1-0.1    0.4396   0.3000    0.0554]); %subplot(514);
% patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X(:,i),...
%     'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
%     'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
% 
% set(ax,'Clim',[mn mx]);
% view([90 0])
% axis(ax,'equal','vis3d','tight');
% camlight(0,180)
% camlight(0,0)
% axis(ax,'on');
% colormap(bipolar(512, 0.99))
% 
% % ax.Position([1 3 4]) = [0.18 0.2350    0.094];
% % ax.Position = [0.1500    0.2396    0.3000    0.0554];
% % ax.Position([3 4]) = [0.75 0.06];
% %%
% % title('front view');
% ax = axes('Position',[ (1/ncols)*i*1-0.1    0.3014    0.3000    0.0554]); %subplot(515);
% patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X(:,i),...
%     'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
%     'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
% set(ax,'Clim',[mn mx]);
% view([-90 0])
% axis(ax,'equal','vis3d','tight');
% camlight(0,180)
% camlight(0,0)
% axis(ax,'on');
% colormap(bipolar(512, 0.99))
% 
% % ax.Position([1 3 4]) = [0.18 0.2350    0.094];
% % ax.Position = [0.1500    0.1014    0.3000    0.0554];
% % ax.Position([3 4]) = [0.75 0.06];
% % title('back view');
% % colorbar
% 
% text(0.05,0.45,'front','units','normalized','Parent', AxesH);
% text(0.05,0.33,'back','units','normalized','Parent', AxesH);
