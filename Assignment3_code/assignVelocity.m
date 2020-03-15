
function [xVelocity,yVelocity,zVelocity] = assignVelocity(varargin)
%assignVelocity Assigns random component velocities given a vector total
%   assignVelocity(system, electron, numDimensions, distribution)
%   Inpus:
%       system - physical system electrons are being simulated in
%       electron - information about electrons being simulated
%       numDimensions - Number of output vectors for 1D, 2D, 3D
%       distribution - type of distribution the velocities will have
%   Outputs:
%       xVelocity  - Vector of x-component velocities
%       yVelocity  - Vector of y-component velocities
%       zVelocity  - Vector of z-component velocities
c.boltzmann = 1.381E-23; %J/K


if nargin == 4 
    
    system = varargin{1};
    electron = varargin{2};
    numDimensions = varargin{3};
    distribution = varargin{4};   
    
    
    %Ignore 1 and 3D cases for now

    %Generate vector of angles, then use angle to find component of x and y

    if(strcmp(distribution, 'normal'))
        %Normally distributed  
    elseif(strcmp(distribution, 'mb'))
        %Maxwell-Boltzmann
        %Using info from here:https://www.mathworks.com/help/matlab/math/random-numbers-with-specific-mean-and-variance.html
        %Assume MB is very close to normal with different std dev and mean
        stddev = sqrt(c.boltzmann.*system.Temp./electron.effM);
        speed = stddev.*randn(1,electron.num) + system.thermalV;

        angleVect = rand(1,electron.num).*2.*pi;
        xVelocity = speed.*cos(angleVect);
        yVelocity = speed.*sin(angleVect);
        zVelocity = [];

        %figure(2)
        %histogram(speed,20);
    else
        angleVect = rand(1,electron.num)*2*pi;
        xVelocity = system.thermalV*cos(angleVect);
        yVelocity = system.thermalV*sin(angleVect);
        zVelocity = [];
    end
elseif nargin == 2
    
    system = varargin{1};
    electron = varargin{2};
        
    stddev = sqrt(c.boltzmann.*system.Temp./electron.effM);
    speed = stddev.*randn() + system.thermalV;

    angleVect = rand().*2.*pi;
    xVelocity = speed.*cos(angleVect);
    yVelocity = speed.*sin(angleVect);
    zVelocity = [];
else
    error("Unknown use case");
end