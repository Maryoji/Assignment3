function [xVect, yVect] = Position(xBound, yBound, numE)
%assignPosition generates random x and y coordinates
%   assignPosition(xBound, yBound, numE)
%   Inputs:
%       xBound - Maximum x value allowed by the system
%       yBound - Maximum y value allowed by the system
%       numE   - Number of electrons being simulated
%   Outputs:
%       xVect  - Vector of x values of length numE
%       yVect  - Vector of y values of length numE

xVect = rand(1,numE)*xBound;
yVect = rand(1,numE)*yBound;

end