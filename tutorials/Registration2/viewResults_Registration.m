addpath '../../utilities/matlab/'

% load source data
% SourcePts contains the coordinates of the vertices, SourceEdges contains a connectivity matrix for the definition of the edges
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

% load deformation of source data
Pts = cell(1,20);
for t=1:20
	Pts{t} = cell(1,6);
	for k=1:6
		Pts{t}{k} = VTKPolyDataReader(['skull_australopithecus_curve_' num2str(k) '__t_' num2str(t-1) '.vtk']);
	end
end

% set color map
color = {'b','r','c','g','y','m'};


% see source and target
figure;
set(gcf,'OuterPosition',[0 0 500 500]);
for t=1:20
	clf;
	hold on
	for p=1:6
		for k=1:size(TargetEdges{p},1)
			plot(TargetPts{p}(TargetEdges{p}(k,:),1),TargetPts{p}(TargetEdges{p}(k,:),2),['-' color{p}],'LineWidth',2);
		end
		for k=1:size(SourceEdges{p},1)
			plot(SourcePts{p}(SourceEdges{p}(k,:),1),SourcePts{p}(SourceEdges{p}(k,:),2),[':' color{p}],'LineWidth',2);
			plot(Pts{t}{p}(SourceEdges{p}(k,:),1),Pts{t}{p}(SourceEdges{p}(k,:),2),['-' color{p}],'LineWidth',2);
		end
	end
	axis([-150 100 -100 120]);
	pause(0.5);
end

figure;
set(gcf,'OuterPosition',[500 0 500 500]);
hold on
for p=1:6
	for k=1:size(TargetEdges{p},1)
		plot(TargetPts{p}(TargetEdges{p}(k,:),1),TargetPts{p}(TargetEdges{p}(k,:),2),'-r','LineWidth',2);
	end
	for k=1:size(SourceEdges{p},1)
		plot(Pts{t}{p}(SourceEdges{p}(k,:),1),Pts{t}{p}(SourceEdges{p}(k,:),2),'-b','LineWidth',2);
	end
end
axis([-150 100 -100 120]);



figure;
set(gcf,'OuterPosition',[0 500 1500 500]);
subplot(1,3,1)
hold on
for p = 1:6
	for k=1:size(SourceEdges{p},1)
		plot(SourcePts{p}(SourceEdges{p}(k,:),1),SourcePts{p}(SourceEdges{p}(k,:),2),'-r','LineWidth',2);
	end
end
axis([-150 100 -100 120]);
subplot(1,3,2)
hold on
for p = 1:6
	for k=1:size(SourceEdges{p},1)
		plot(Pts{20}{p}(SourceEdges{p}(k,:),1),Pts{20}{p}(SourceEdges{p}(k,:),2),'-b','LineWidth',2);	
	end
end
axis([-150 100 -100 120]);
subplot(1,3,3)
hold on
for p = 1:6
	for k=1:size(TargetEdges{p},1)
		plot(TargetPts{p}(TargetEdges{p}(k,:),1),TargetPts{p}(TargetEdges{p}(k,:),2),'-r','LineWidth',2);
	end
end
axis([-150 100 -100 120]);



