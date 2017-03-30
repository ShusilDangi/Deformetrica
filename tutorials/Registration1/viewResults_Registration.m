addpath '../../utilities/matlab/'

% load source data
% SourcePts contains the coordinates of the vertices, SourceEdges contains a connectivity matrix for the definition of the edges
[SourcePts SourceEdges] = VTKPolyDataReader(['skull_australopithecus.vtk']);
% load target data
[TargetPts TargetEdges] = VTKPolyDataReader(['../skull_sapiens.vtk']);

% load deformation of source data
for t=1:20
	Pts{t} = VTKPolyDataReader(['skull_australopithecus__t_' num2str(t-1) '.vtk']);
end


% see source and target
figure;
set(gcf,'OuterPosition',[0 0 500 500]);
for t=1:20
	clf;
	hold on
	for k=1:size(TargetEdges,1)
		plot(TargetPts(TargetEdges(k,:),1),TargetPts(TargetEdges(k,:),2),'-r','LineWidth',2);
	end
	for k=1:size(SourceEdges,1)
		plot(SourcePts(SourceEdges(k,:),1),SourcePts(SourceEdges(k,:),2),'-r','LineWidth',2);
		plot(Pts{t}(SourceEdges(k,:),1),Pts{t}(SourceEdges(k,:),2),'-b','LineWidth',2);
	end
	axis([-150 100 -100 120]);
	pause(0.5);
end


figure;
set(gcf,'OuterPosition',[0 500 1500 500]);
subplot(1,3,1)
hold on
for k=1:size(SourceEdges,1)
	plot(SourcePts(SourceEdges(k,:),1),SourcePts(SourceEdges(k,:),2),'-r','LineWidth',2);
end
axis([-150 100 -100 120]);
subplot(1,3,2)
hold on
for k=1:size(SourceEdges,1)
	plot(Pts{20}(SourceEdges(k,:),1),Pts{20}(SourceEdges(k,:),2),'-b','LineWidth',2);	
end
axis([-150 100 -100 120]);
subplot(1,3,3)
hold on
for k=1:size(TargetEdges,1)
	plot(TargetPts(TargetEdges(k,:),1),TargetPts(TargetEdges(k,:),2),'-r','LineWidth',2);
end
axis([-150 100 -100 120]);



