addpath '../../utilities/matlab/'

% SourcePts contains the coordinates of the vertices, SourceEdges contains a connectivity matrix for the definition of the edges
SourcePts = cell(1,6);
SourceEdges = cell(1,6);
for k=1:6
	[SourcePts{k} SourceEdges{k}] = VTKPolyDataReader(['skull_australopithecus_curve_' num2str(k) '.vtk']);
end

% load deformation of source data
Pts = cell(1,20);
Lat = cell(1,20);
for t=1:20
	Pts{t} = cell(1,6);
	for k=1:6
		Pts{t}{k} = VTKPolyDataReader(['skull_australopithecus_curve_' num2str(k) '__t_' num2str(t-1) '.vtk']);
	end
	[Lat{t} LatticeEdges] = VTKPolyDataReader(['lattice_flow__t_' num2str(t-1) '.vtk']);
end

% set color map
color = {'b','r','c','g','y','m'};

figure;
set(gcf,'OuterPosition',[0 0 500 500]);
for t=1:20
	clf;
	hold on
	for k=1:size(LatticeEdges,1)
		plot(Lat{t}(LatticeEdges(k,:),1),Lat{t}(LatticeEdges(k,:),2),'-k','LineWidth',2);
	end
	for p=1:6
		for k=1:size(SourceEdges{p},1)
			plot(Pts{t}{p}(SourceEdges{p}(k,:),1),Pts{t}{p}(SourceEdges{p}(k,:),2),color{p},'LineWidth',2);
		end
	end
	axis([-150 100 -100 120]);
	pause(0.5);
end
