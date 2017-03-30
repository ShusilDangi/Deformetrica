addpath '../../utilities/matlab/'

% load source data
% SourcePts contains the coordinates of the vertices, SourceEdges contains a connectivity matrix for the definition of the edges
[AP AE] = VTKPolyDataReader('skull_australopithecus.vtk');
[NP NE] = VTKPolyDataReader('skull_neandertalis.vtk');

% load deformation of source data and lattice deformation
S = cell(1,20);
N = cell(1,20);
for t=1:20
	[S{t} E] = VTKPolyDataReader(['deformed_australopithecus_flow__t_' num2str(t-1) '.vtk']);
	N{t} = VTKPolyDataReader(['skull_neandertalis_flow__t_' num2str(t-1) '.vtk']);
end

% see source and target
figure;
set(gcf,'OuterPosition',[0 0 500 500]);
for t=1:20
	clf;
	hold on
	for k=1:size(AE,1)
		plot(AP(AE(k,:),1),AP(AE(k,:),2),'-b','LineWidth',2);
	end
	for k=1:size(E,1)
		plot(S{t}(E(k,:),1),S{t}(E(k,:),2),'-r','LineWidth',2);
	end
	axis([-150 100 -100 120]);
	pause(0.5);
end



figure;
set(gcf,'OuterPosition',[500 0 500 500]);
for t=1:20
	clf;
	hold on
	for k=1:size(AE,1)
		plot(AP(AE(k,:),1),AP(AE(k,:),2),'-b','LineWidth',2);
	end
	for k=1:size(NE,1)
		plot(N{t}(NE(k,:),1),N{t}(NE(k,:),2),'-r','LineWidth',2);
	end
	axis([-150 100 -100 120]);
	pause(0.5);
end
