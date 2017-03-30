addpath '../../utilities/matlab/'

% load source data
% SourcePts contains the coordinates of the vertices, SourceEdges contains a connectivity matrix for the definition of the edges
[SourcePts SourceEdges] = VTKPolyDataReader('skull_australopithecus.vtk');
[LatticePts LatticeEdges] = VTKPolyDataReader('lattice.vtk');

% load deformation of source data and lattice deformation
Pts = cell(1,20);
Lat = cell(1,20);
for t=1:20
	Pts{t} = VTKPolyDataReader(['skull_australopithecus__t_' num2str(t-1) '.vtk']);
	Lat{t} = VTKPolyDataReader(['lattice_flow__t_' num2str(t-1) '.vtk']);
end


% see source and target
figure;
set(gcf,'OuterPosition',[0 0 500 500]);
for t=1:20
	clf;
	hold on
	for k=1:size(LatticeEdges,1)
		plot(Lat{t}(LatticeEdges(k,:),1),Lat{t}(LatticeEdges(k,:),2),'-k','LineWidth',2);
	end
	for k=1:size(SourceEdges,1)
		plot(Pts{t}(SourceEdges(k,:),1),Pts{t}(SourceEdges(k,:),2),'-b','LineWidth',2);
	end
	axis([-150 100 -100 120]);
	pause(0.5);
end
