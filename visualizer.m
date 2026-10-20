%% Clear stuff
clear; clc;
close all;

%% Info
itt = 400; % number of itterations
L = 250; % length of domain
W = L; % width of domain

framerate = 1; % fps of videos 
videoDimensions = [300, 300, 652, 514]; % dimensions of videos

%% Load Data
%rhoCpp = loadData("cpp/output/rho", itt);
u = readmatrix("cpp/output/velocity/uX.dat");
rho = readmatrix("cpp/output/density/rhoX.dat");

ghia = load("Ghia_Re100.mat");

ux = u(1:W, 1:L);
uy = u(W+1:W*2, 1:L);
%%

figure;
hold on;
plot((1:W)/W, ux(:,ceil(W/2))/max(abs(ux(:,ceil(W/2)))), "y");
plot(ghia.y_Ghia, ghia.u_Ghia/(max(abs(ghia.u_Ghia))), "--y");
plot((1:L)/L, uy(ceil(L/2),:)/max(abs(uy(ceil(L/2),:))), "b");
plot(ghia.x_Ghia, ghia.v_Ghia/max(abs(ghia.v_Ghia)), "--b");
axis equal tight;
box on;
hold off;

figure;
hold on;
contourf(rho);
axis equal tight;
box on;
hold off;

figure;
hold on;
quiver(ux,uy);
axis equal tight;
box on;
hold off;
%% Export fluid density
%minmaxRhoCpp = [getMin(rhoCpp, itt) getMax(rhoCpp, itt)];

%%
%MP4ExportContourf("rhoCpp", framerate, videoDimensions, rhoCpp, itt, NaN, {1:L, 1:W}, minmaxRhoCpp, 30);


%% Load function
function [out] = loadData(path, itt)

    out = cell(1, itt+1);
        out{1} = readmatrix(path+"0.dat");

    for i = 1:itt
        filename = path+(i)+".dat";
        out{i+1} = readmatrix(filename);
        % if mod(i, 100) == 0
        %     fprintf("Table %i read \n", i);
        % end
    end
   
    fprintf("Imported data from: " + path + "\n");
end

%% MP4 Export
function MP4ExportMesh(vidTitle, data, L, W, itt, dt, limFlag, lim)
    % Initialize video writer
    newVid = VideoWriter("videos/"+vidTitle+".mp4", 'MPEG-4');
    newVid.FrameRate = 30;
    open(newVid);
    
    % Create figure once
    fig = figure('Position', [100, 100, 652, 514]);  % Set consistent size
    axis tight manual
    
    %guh = [min(data{3},[], "all"), max(data{1},[], "all")];

    % Loop through and render each mesh
    fprintf("loop started");
    for i = 1:(itt+1)
        Z = data{i};
        mesh(1:L, 1:W, Z);

        if limFlag
            zlim(lim)
        end
   
        title(sprintf('Frame %d; %.3e sec', i, dt*(i-1)));
        drawnow;
     
        % Capture frame from current figure
        frame = getframe(fig);
        writeVideo(newVid, frame);
        if mod(i, 100) == 0
            fprintf("Table %i wrote \n", i);
        end
    end
    
    close(newVid);
    close(fig);
end

%% MP4 Export
function MP4ExportContourf(vidTitle, framerate, vidDim, data, itt, dt, plotDim, minmax, detail)
    % Initialize video writer
    newVid = VideoWriter(vidTitle+".mp4", 'MPEG-4');
    newVid.FrameRate = framerate;
    open(newVid);
    
    % Create figure once
    fig = figure('Position', vidDim);  % Set consistent size
    axis tight manual

    % Loop through and render each mesh
    fprintf("loop started");
    for i = 1:(itt+1)
        Z = data{i};
        [M, c] = contourf(Z(plotDim{1}, plotDim{2}),'LineColor','none');
        c.LevelListMode = "manual";
        c.LevelList = linspace(minmax(1), minmax(2), detail);
        colorbar;

        title(sprintf('Frame %d; %.3e sec', i, dt*(i-1)));

        drawnow;
     
        % Capture frame from current figure
        frame = getframe(fig);
        writeVideo(newVid, frame);
    end
    
    % Close video and figure objects
    close(newVid);
    close(fig);

    fprintf(vidTitle+".mp4 exported\n");
end

%% Process velocity function
function [out] = uMag(u, L, itt)

    out = cell(1, itt+1);

    for i = 1:itt+1
        X = u{i};
        out{i} = sqrt(X(1:L, :).*X(1:L, :) + X(L+1:2*L, :).*X(L+1:2*L, :));
    end

end

%% Min
function [tMin] = getMin(A, itt)

    tMin = min(A{1},[],"all");

    for i = 2:itt+1
        cMin = min(A{i},[],"all");
        if tMin > cMin
            tMin = cMin;
        end
    end

end

%% Max
function [tMax] = getMax(A, itt)

    tMax = max(A{1},[],"all");

    for i = 2:itt+1
        cMax = max(A{i},[],"all");
        if tMax < cMax
            tMax = cMax;
        end
    end

end