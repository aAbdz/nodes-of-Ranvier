%% run_nodes_of_ranvier_pipeline.m

% Requirements:
%   - Bio-Formats MATLAB toolbox
%   - bfOpen3DVolume
%   - util_forcept_segmentation
%   - util_plot
%   - util_nr_quantification
%   - optional: imshow3D

clear; clc;

%% ---------------- User settings ----------------

cfg.codeDir = '';
cfg.inputDir = '';
cfg.outputDir = '';

cfg.voxelSize = [0.04, 0.04, 0.1];   % [x y z] in micrometers
cfg.referenceXY = 0.04;              % used for anisotropy correction
cfg.medianFilterSize = [5, 5, 5];

cfg.displayResults = false;
cfg.saveFigures = true;
cfg.overwrite = true;

%% ---------------- Run pipeline ----------------

addpath(genpath(cfg.codeDir));

if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end

sampleDirs = list_subdirectories(cfg.inputDir);

for iSample = 1:numel(sampleDirs)

    sampleName = sampleDirs(iSample).name;
    sampleInputDir = fullfile(cfg.inputDir, sampleName);
    sampleOutputDir = fullfile(cfg.outputDir, sampleName);

    if ~exist(sampleOutputDir, 'dir')
        mkdir(sampleOutputDir);
    end

    imageFiles = list_image_files(sampleInputDir);

    fprintf('\nProcessing sample %d/%d: %s\n', ...
        iSample, numel(sampleDirs), sampleName);

    for iFile = 1:numel(imageFiles)

        fileName = imageFiles(iFile).name;
        [~, baseName, ~] = fileparts(fileName);

        fprintf('  Processing file %d/%d: %s\n', ...
            iFile, numel(imageFiles), fileName);

        try
            process_single_volume(sampleInputDir, sampleOutputDir, fileName, baseName, cfg);
        catch ME
            warning('Failed to process %s: %s', fileName, ME.message);
        end
    end
end

fprintf('\nDone.\n');

%% ========================================================================
%% Local functions
%% ========================================================================

function process_single_volume(inputDir, outputDir, fileName, baseName, cfg)

    outputMatFile = fullfile(outputDir, [baseName, '.mat']);

    if exist(outputMatFile, 'file') && ~cfg.overwrite
        fprintf('    Skipping existing result: %s\n', outputMatFile);
        return;
    end

    imagePath = fullfile(inputDir, fileName);

    % Read original volume
    rawVolumeCell = bfOpen3DVolume(imagePath);
    rawVolume = rawVolumeCell{1}{1};

    % Find matching Airyscan-processed file
    airyscanFile = find_airyscan_file(inputDir, baseName);

    if isempty(airyscanFile)
        fprintf('    No matching Airyscan file found. Skipping.\n');
        return;
    end

    airyscanPath = fullfile(inputDir, airyscanFile);
    airyscanCell = bfOpen3DVolume(airyscanPath);
    airyscanVolume = airyscanCell{1}{1};

    % Infer number of z-slices from Airyscan file
    nZ = size(airyscanVolume, 3) / 4;
    nZ = floor(nZ);

    % Extract green channel from interleaved channels
    greenVolume = rawVolume(:, :, 2:4:end);
    greenVolume = greenVolume(:, :, 1:nZ);

    % Normalize to uint8
    greenVolume = normalize_to_uint8(greenVolume);

    % Enhance image
    filteredVolume = medfilt3(greenVolume, cfg.medianFilterSize);
    maxProjection = max(filteredVolume, [], 3);

    % Segment nodes of Ranvier
    voxelRatio = cfg.voxelSize / cfg.referenceXY;
    labelVolume = util_forcept_segmentation(maxProjection, filteredVolume, voxelRatio);

    % Display optional quality control
    if cfg.displayResults
        show_segmentation(labelVolume);
    end

    % Save figure
    if cfg.saveFigures
        figureBaseName = fullfile(outputDir, baseName);
        util_plot(maxProjection, labelVolume, figureBaseName);
    end

    % Quantification
    quants_pn = util_nr_quantification(labelVolume, cfg.voxelSize);

    % Save result
    lbl = labelVolume; %#ok<NASGU>
    save(outputMatFile, 'lbl', 'quants_pn', 'cfg', '-v7.3');

    close all;
end

function subdirs = list_subdirectories(parentDir)

    allItems = dir(parentDir);
    isValid = [allItems.isdir] & ~ismember({allItems.name}, {'.', '..'});
    subdirs = allItems(isValid);
end

function files = list_image_files(folderPath)

    allItems = dir(folderPath);
    allItems = allItems(~[allItems.isdir]);

    names = {allItems.name};

    skipMask = contains(names, 'Processing', 'IgnoreCase', true) | ...
               contains(names, 'Thumbs.db', 'IgnoreCase', true) | ...
               contains(names, 'Airyscan', 'IgnoreCase', true);

    files = allItems(~skipMask);
end

function airyscanFile = find_airyscan_file(folderPath, baseName)

    safeBaseName = strrep(baseName, '.', '_');
    searchPattern = fullfile(folderPath, [safeBaseName, '-Airyscan*']);

    matches = dir(searchPattern);

    if isempty(matches)
        airyscanFile = '';
    else
        airyscanFile = matches(1).name;
    end
end

function out = normalize_to_uint8(volume)

    volume = double(volume);
    maxValue = max(volume(:));

    if maxValue == 0
        out = uint8(volume);
    else
        out = uint8(255 * volume / maxValue);
    end
end

function show_segmentation(labelVolume)
    figure;
    p = patch(isosurface(smooth3(logical(labelVolume), 'gaussian', 7)));
    hold on;

    p.FaceColor = [0.1, 0.2, 0.1];
    p.FaceAlpha = 0.05;
    p.EdgeColor = [0.1, 0.2, 0.1];
    p.EdgeAlpha = 0.3;

    camlight;
    axis equal;
    set(gca, 'Visible', 'off', 'Color', 'k');
    drawnow;
end