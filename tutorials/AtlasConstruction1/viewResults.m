addpath '../../utilities/matlab/'

% load estimated template shape
[TemplatePts TemplateEdges] = VTKPolyDataReader(['skull_template.vtk']);

% display final template
figure;
hold on
for k=1:size(TemplateEdges,1)
	plot(TemplatePts(TemplateEdges(k,:),1),TemplatePts(TemplateEdges(k,:),2),'-b','LineWidth',3);
end
axis([-150 100 -100 120]);

% load optimal position of control points
CP = load('CP_final.txt');

% load set of optimal momentum vectors
MOM = readMomentaFile('MOM_final.txt');

name = {'australopithecus','habilis','erectus','neandertalis','sapiens'};

% load data
P = cell(1,5);
E = cell(1,5);
for s = 1:5
	[P{s} E{s}] = VTKPolyDataReader(['../skull_' name{s} '.vtk']);
end

% load the frames of the 5 template-to-subjects deformations
Pt = cell(1,20);
for t=1:20
	Pt{t} = cell(1,5);
	for s=1:5
		Pt{t}{s} = VTKPolyDataReader(['skull_to_subject_' num2str(s-1) '__t_' num2str(t-1) '.vtk']);
	end
end



% display template-to-subject deformations
color = {'g','c','m','k','y'};
pos = [3 4 1 2 6];
figure;
set(gcf,'OuterPosition',[0 0 1500 1000]);
for t=1:20
	clf;
	for s=1:5
		subplot(2,3,pos(s));
		hold on
		for k=1:size(E{s},1)
			plot(P{s}(E{s}(k,:),1),P{s}(E{s}(k,:),2),'-r','LineWidth',3);
		end
		for k=1:size(TemplateEdges,1)
			plot(Pt{t}{s}(TemplateEdges(k,:),1),Pt{t}{s}(TemplateEdges(k,:),2),color{s},'LineWidth',3);
		end
		for k=1:size(TemplateEdges,1)
			plot(TemplatePts(TemplateEdges(k,:),1),TemplatePts(TemplateEdges(k,:),2),'-b','LineWidth',3*(t==1)+(t>1));
		end
		title(name{s});
		axis([-150 100 -100 120]);
	end
	subplot(2,3,5)
	hold on
	for s=1:5
		quiver(CP(:,1),CP(:,2),MOM(1,:,s)',MOM(2,:,s)',0.5,color{s},'LineWidth',3);
	end
	for k=1:size(TemplateEdges,1)
		plot(TemplatePts(TemplateEdges(k,:),1),TemplatePts(TemplateEdges(k,:),2),'-b','LineWidth',3);
	end
	title('Estimated template and momentum vectors');
	axis([-150 100 -100 120]);
	pause(0.5);
end
