function part3(system, electron, numLoops, timeStep, boundType, var1, var2)
%noCollision Simulates electrons moving through silicon with collision and
%a bottleneck region
%   bottleNeck(system, electron, numLoops, timeStep)
%   Inpus:
%       system     - Structure containing properties of silicon sample
%       electron   - Structure containing properties of electrons
%       numLoops   - Number of loops that the simulation will run for
%       timeStep   - Length of time between each plot update
%       boundType  - String denoting barrier heaviour, either 'specular' or
%       'diffusive'
%       varargin1  - X component of electric field
%       varargin2  - Y component of electric field 
%   Outputs:
%       None

Ex = var1;
Ey = var2;

%Define forbidden region 'boxes'
boxydim = [0.4E-7, 0.6E-7];
boxxdim = [0.8E-7, 1.2E-7];

%Draw box 1
figure('Renderer', 'painters', 'Position', [10 10 1100 600])
hold on
plot([boxxdim(1) boxxdim(1)], [system.y boxydim(2)], '-k');
plot([boxxdim(1) boxxdim(2)], [boxydim(2) boxydim(2)], '-k');
plot([boxxdim(2) boxxdim(2)], [system.y boxydim(2)], '-k');

%Draw box 2
plot([boxxdim(1) boxxdim(1)], [0 boxydim(1)], '-k');
plot([boxxdim(1) boxxdim(2)], [boxydim(1) boxydim(1)], '-k');
plot([boxxdim(2) boxxdim(2)], [0 boxydim(1)], '-k');

%Move electrons that start existing in the boxes
for i = 1:electron.num
    while((electron.x(1,i) >= boxxdim(1) && electron.x(1,i) <= boxxdim(2)) && (electron.y(1,i) <= boxydim(1) || electron.y(1,i) >= boxydim(2)))
        [electron.x(1,i),electron.y(1,i)] = Position(system.x, system.y, 1);
    end
end



xlim([0 system.x]);
ylim([0 system.y]);

c.boltzmann = 1.381E-23; %J/K
t = [0];

%Scatter probability from assignment outline
pScatter = 1 - exp(-timeStep/system.Tmn);

[electron.vx, electron.vy] = assignVelocity(system, electron, 2, 'uniform');

colourScheme = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];


for i = 1:numLoops
    %Assuming x and y are the same size, find the length for boundary
    %conditions
    len = size(electron.x, 1);
    t(len) = timeStep.*i;
    
    jumped(len,:) = zeros(1,electron.num);
    
    %Boundary conditions and scattering calcs
    for j = 1:electron.num        
        if(strcmp('specular', boundType))
            %Specular boundary condition case
            
            %First handle boundaries
            if(electron.x(len, j) < 0)
                electron.vx(len,j) = -electron.vx(len,j);
                electron.x(len, j) = 0;
            elseif(electron.x(len,j) > system.x)
                electron.vx(len,j) = -electron.vx(len,j);
                electron.x(len, j) = system.x;
            end
                
            if(electron.y(len,j) < 0)
                electron.vy(len,j) = -electron.vy(len,j);
                electron.y(len, j) = 0;
            elseif(electron.y(len,j) > system.y)
                electron.vy(len,j) = -electron.vy(len,j);
                electron.y(len, j) = system.y;
            end
                
            
            if electron.y(len,j) < boxydim(1) %Check box 1
                if (electron.x(len,j) >= boxxdim(1) && electron.x(len,j) <= boxxdim(2))
                    %Electron is inside top square
                    %Check two cycles ago to find where the electron came
                    %from
                    if(electron.x(len-1,j) <= boxxdim(1) || electron.x(len-1,j) >= boxxdim(2))
                        electron.vx(len,j) = -electron.vx(len,j);
                        electron.x(len,j) = electron.x(len-1,j);
                        electron.y(len,j) = electron.y(len-1,j);
                    else
                        electron.vy(len,j) = -electron.vy(len,j);
                        electron.x(len,j) = electron.x(len-1,j);
                        electron.y(len,j) = electron.y(len-1,j);
                    end 
                end
            elseif electron.y(len,j) > boxydim(2) %Check box 2
                if (electron.x(len,j) >= boxxdim(1) && electron.x(len,j) <= boxxdim(2))
                    %Electron is inside bottom square
                    if(electron.x(len-1,j) <= boxxdim(1) || electron.x(len-1,j) >= boxxdim(2))
                        electron.vx(len,j) = -electron.vx(len,j);
                        electron.x(len,j) = electron.x(len-1,j);
                        electron.y(len,j) = electron.y(len-1,j);
                    else
                        electron.vy(len,j) = -electron.vy(len,j);
                        electron.x(len,j) = electron.x(len-1,j);
                        electron.y(len,j) = electron.y(len-1,j);
                    end 
                end
            end
            

            
            
        elseif(strcmp('diffusive', boundType))
            %diffusive (rethermalization) case
            
            %First handle boundaries
            
            %First handle boundaries
            if(electron.x(len,j) < 0)
                [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
                electron.x(len, j) = 0;
            elseif(electron.x(len,j) > system.x)
                [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
                electron.x(len, j) = system.x;
            end
                
            if(electron.y(len,j) < 0)
                [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
                electron.y(len, j) = 0;
            elseif(electron.y(len,j) > system.y)
                [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
                electron.y(len, j) = system.y;
            end
                
            
            if electron.y(len,j) < boxydim(1) %Check box 1
                if (electron.x(len,j) >= boxxdim(1) && electron.x(len,j) <= boxxdim(2))
                    %Electron is inside top square
                    %Check two cycles ago to find where the electron came
                    %from
                    if(electron.x(len-1,j) <= boxxdim(1) || electron.x(len-1,j) >= boxxdim(2))
                        [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
                        electron.x(len,j) = electron.x(len-1,j);
                        electron.y(len,j) = electron.y(len-1,j);
                    else
                        [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
                        electron.x(len,j) = electron.x(len-1,j);
                        electron.y(len,j) = electron.y(len-1,j);
                    end 
                end
            elseif electron.y(len,j) > boxydim(2) %Check box 2
                if (electron.x(len,j) >= boxxdim(1) && electron.x(len,j) <= boxxdim(2))
                    %Electron is inside bottom square
                    if(electron.x(len-1,j) <= boxxdim(1) || electron.x(len-1,j) >= boxxdim(2))
                        [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
                        electron.x(len,j) = electron.x(len-1,j);
                        electron.y(len,j) = electron.y(len-1,j);
                    else
                        [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
                        electron.x(len,j) = electron.x(len-1,j);
                        electron.y(len,j) = electron.y(len-1,j);
                    end 
                end
            end
            
        else
            error('Unknown boundary type')
        end
        
        
        %pScatter is quite small, checking for less than pScatter gives about
        %3% proability
        if(rand() < pScatter)         
            %Rethermalize using overladed assignVelocity
            [electron.vx(len,j), electron.vy(len,j)] = assignVelocity(system, electron);
            
        end      
    end
    
    
    
    %Plot all of the electrons if there are fewer than 7
    %Uncomment to plot during loop 
    % for j = 1:min([electron.num 7])
    %     noLine = find(jumped(:,j),50);
    %     if isempty(noLine)
    %         %No skips, plot normally
    %         plot(electron.x(:,j),electron.y(:,j), 'Color', colourScheme(j,:));
    %     else
    %         plot(electron.x((1:noLine(1)-1),j),electron.y((1:noLine(1)-1),j), 'Color', colourScheme(j,:));
    %         for k = 1:length(noLine)         
    %             if k == length(noLine)
    %                 plot(electron.x((noLine(k):end),j),electron.y((noLine(k):end),j), 'Color', colourScheme(j,:));
    %             else
    %                 plot(electron.x((noLine(k):noLine(k+1)-1),j),electron.y((noLine(k):noLine(k+1)-1),j), 'Color', colourScheme(j,:));
    %             end
    %         end
    %     end
    % end
    
    % pause(0.01);
    
    % [electron.x, electron.vx] = propagate(electron.x, electron.vx, timeStep);
    % [electron.y, electron.vy] = propagate(electron.y, electron.vy, timeStep);
    electron = Propagate(electron, system, timeStep, Ex, Ey);

    
 %Plotting after loop, faster
end
 

 
 for j = 1:min([electron.num 7])
    noLine = find(jumped(:,j),50);
    if isempty(noLine)
        plot(electron.x(:,j),electron.y(:,j), 'Color', colourScheme(j,:));
    else
        plot(electron.x((1:noLine(1)-1),j),electron.y((1:noLine(1)-1),j), 'Color', colourScheme(j,:));
        for k = 1:length(noLine)         
            if k == length(noLine)
                plot(electron.x((noLine(k):end),j),electron.y((noLine(k):end),j), 'Color', colourScheme(j,:));
            else
                plot(electron.x((noLine(k):noLine(k+1)-1),j),electron.y((noLine(k):noLine(k+1)-1),j), 'Color', colourScheme(j,:));
            end
        end
    end
 end 

title("Electron Trajectories With Calculated Electric Field", 'Interpreter', 'Latex')
xlabel("X Dimension (nm)", 'Interpreter', 'Latex');
ylabel("Y Dimension (nm)", 'Interpreter', 'Latex');
set(gca, 'FontSize', 12);

 
 
%Map generation

%Density, temperature maps
xdim = linspace(0, system.x, 30);
ydim = linspace(0,system.y, 15);
zDense = zeros(length(xdim),length(ydim));

xTempBin = discretize(electron.x(end,:), xdim);
yTempBin = discretize(electron.y(end,:), ydim);
figure('Renderer', 'painters', 'Position', [10 10 1100 600])
for i = 1:length(xdim)
    for j = 1:length(ydim)
        %Translate logical to indices of electrons in region
        inRegion = (xTempBin == i) & (yTempBin == j);
        inRegion = find(inRegion, electron.num);    
        
        %Count number of electrons for density maps
        zDense(i,j) = length(inRegion);        
    end
end

%Density plot
[X,Y] = meshgrid(xdim, ydim);
figure('Renderer', 'painters', 'Position', [10 10 1100 600])
surf(X,Y, zDense');
title("Electron Density", 'interpreter', 'latex');
xlabel('X Dimension (m)', 'interpreter', 'latex');
ylabel('Y Dimension (m)', 'interpreter', 'latex');
set(gca, 'fontsize', 15)
hold on
    
    




end
