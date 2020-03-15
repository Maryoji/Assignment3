function [electron] = Propagate(electron, system, timeStep, varargin1, varargin2)
%propagate Returns a new electron object with updated positions based on Newton's laws
%   propagate(location, velocity, system, timeStep)
%   Inpus:
%       location - vector marking original locations
%       velocity - vector marking speed of particles
%       timeStep - Amount of time particle is travelling before update
%   Outputs:
%       newLocation - Vector containing new locations

q = 1.602E-19; %C

Ex = varargin1;
Ey = varargin2;

%Assumed vx and vy have same length
currentSize = size(electron.vx, 1);

%Add a row of blank zeros to get correct size for new locations
electron.x(currentSize + 1,:) = zeros(1,size(electron.x,2));
electron.vx(currentSize + 1,:) = zeros(1,size(electron.x,2));

electron.y(currentSize + 1,:) = zeros(1,size(electron.y,2));
electron.vy(currentSize + 1,:) = zeros(1,size(electron.y,2));

%Calculate the acceleration due to the presence of an efield

if nargin == 3
	fieldForceX = q.*system.EfieldX;
	fieldAccelX = fieldForceX./electron.effM;

	fieldForceY = q.*system.EfieldY;
	fieldAccelY = fieldForceY./electron.effM;
else
    coordinateX = round(electron.x(currentSize,:).*1E9);
    coordinateY = round(electron.y(currentSize,:).*1E9);

    for i = 1:length(coordinateX)
        if coordinateX(i) == 0
            coordinateX(i) = 1;
        end
        if coordinateY(i) == 0
            coordinateY(i) = 1;
        end
        try
            fieldForceX(i) = q.*Ex(coordinateY(i),coordinateX(i));
            fieldForceY(i) = q.*Ey(coordinateY(i),coordinateX(i));	
        catch e
            fprintf(e.message);
        end
    end
    
    fieldAccelX = fieldForceX./electron.effM;
    fieldAccelY = fieldForceY./electron.effM;
end

 


%Add the new location values into the new row, using old currentSize
try    
    electron.x(currentSize + 1,:) = electron.x(currentSize,:) + electron.vx(currentSize,:)*timeStep;
    electron.vx(currentSize + 1,:) = electron.vx(currentSize,:) + fieldAccelX*timeStep;

    electron.y(currentSize + 1,:) = electron.y(currentSize,:) + electron.vy(currentSize,:)*timeStep;
    electron.vy(currentSize + 1,:) = electron.vy(currentSize,:) + fieldAccelY*timeStep;
catch e
    fprintf(e.message);
end

%newLocation = location;

end