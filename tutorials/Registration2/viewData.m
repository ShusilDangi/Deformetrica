addpath '../../utilities/matlab/'

% load source data
% SourcePts contains the coordinates of the vertices, SourceEdges contains a connectivity matrix for the definition of the edges
% source data consists in 6 labelled parts
SourcePts = cell(1,6);
SourceEdges = cell(1,6);
for k=1:6
	[SourcePts{k} SourceEdges{k}] = VTKPolyDataReader(['skull_australopithecus_curve_' num2str(k) '.vtk']);
end



% load target data
TargetPts = cell(1,6);
TargetEdges = cell(1,6);
for k=1:6
	[TargetPts{k} TargetEdges{k}] = VTKPolyDataReader(['skull_sapiens_curve_' num2str(k) '.vtk']);
end

% set color map
color = {'b','r','c','g','y','m'};


% see source and target
figure;
set(gcf,'OuterPosition',[0 0 1000 500]);
subplot(1,2,1)
hold on
for p=1:6
	for k=1:size(SourceEdges{p},1)
		plot(SourcePts{p}(SourceEdges{p}(k,:),1),SourcePts{p}(SourceEdges{p}(k,:),2),[':' color{p}],'LineWidth',2);
	end
end
axis([-150 100 -100 120]);
subplot(1,2,2)
hold on
for p=1:6
	for k=1:size(TargetEdges{p},1)
		plot(TargetPts{p}(TargetEdges{p}(k,:),1),TargetPts{p}(TargetEdges{p}(k,:),2),['-' color{p}],'LineWidth',2);
	end
end
axis([-150 100 -100 120]);
