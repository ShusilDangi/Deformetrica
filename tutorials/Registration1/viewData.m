addpath '../../utilities/matlab/'

% load source data
% SourcePts contains the coordinates of the vertices, SourceEdges contains a connectivity matrix for the definition of the edges
[SourcePts SourceEdges] = VTKPolyDataReader(['../skull_australopithecus.vtk']);
% load target data
[TargetPts TargetEdges] = VTKPolyDataReader(['../skull_sapiens.vtk']);


% see source and target
figure;
subplot(1,2,1)
hold on
for k=1:size(SourceEdges,1)
	plot(SourcePts(SourceEdges(k,:),1),SourcePts(SourceEdges(k,:),2),'-r','LineWidth',3);
end
axis([-150 100 -100 120]);
subplot(1,2,2)
hold on
for k=1:size(TargetEdges,1)
	plot(TargetPts(TargetEdges(k,:),1),TargetPts(TargetEdges(k,:),2),'-r','LineWidth',3);
end
axis([-150 100 -100 120]);
set(gcf,'OuterPosition',[0 0 1000 500]);
