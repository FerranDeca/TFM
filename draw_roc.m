function draw_roc(roc,varargin)
%DRAW_ROC Function used to plot the roc curve
%   roc --> two column matrix containing the curve points
%   first (varargin{1}) --> draw the 45 degree line if 1;
%   type (varargin{2}) --> draw the normal roc if empty or 'roc' or draws
%   the log_roc if 'log'

if length(varargin) <= 1 | varargin{2} == 'roc'

    if  isempty(varargin) | varargin{1} == 1
        P1=[0 0];P2=[1 1];
        line = plot([P1(1) P2(1)],[P1(2) P2(2)],'--');
        set(get(get(line,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off');
        hold on
    end
    plot(roc(:,1),roc(:,2));
    
    
else
    loglog(roc(:,1),1-roc(:,2));
end
ylabel('P_{D}');
xlabel('P_{FA}');
    
end

