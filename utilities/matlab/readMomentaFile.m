%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      %
%                                     Deformetrica                                     %
%                                                                                      %
%    Copyright Inria and the University of Utah.  All rights reserved. This file is    %
%    distributed under the terms of the Inria Non-Commercial License Agreement.        %
%                                                                                      %
%                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MOM = readMomentaFile(filename)

% usage MOM = readMomentaFile(filename)

fid = fopen(filename, 'r');
if (fid==-1)
    fprintf(1,'Error: file descriptor not valid, check the file name.\n');
    return;
end

newL = fgetl(fid);

A = sscanf(newL,'%d %d %d');
nsub = A(1);
nrow = A(2);
ncol = A(3);

MOM = zeros(ncol,nrow,nsub);

for su = 1:nsub
	
	newL = fgetl(fid);
	
	for r = 1:nrow
		newL = fgetl(fid);
		[num c] = sscanf(newL,'%f');
		MOM(:,r,su) = num;
	end

end

fclose(fid);

end
