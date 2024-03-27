% QuaPOS-TEM: QUAntification of mouse Photoreceptor Outer Segments in 
%             Transmission Electron Microscopy images
%
%  Analyse orientations and mutual alignment of membrane stacks in the 
%  outer segments of photoreceptor cells of mouse retina
%
% Karl Hoffmann, Max Planck Institute of Molecular Cell Biology and
% Genetics, Dresden, Germany
% generated: 2020-01-10
% last edit: 2024-03-06
% code last run under MATLAB R2023b

%% folder structure for batch processing
% assuming that head_input_folder contains subfolders named by specimen age, 
% for example 10, 12, 16 and 20 weeks of age.
% timePoints selects subfolders to be analysed. They shall contain images of 
% the same file format ('.bmp' for the publication), each accompagnied with 
% a .txt file for metadata of imaging conditions like resolution, and 
% accompagnied with information on Regions Of Interest (ROIs) stored in a file 
% of essentially the same name (up to leading or trailing "ROI", "ROIs" or 
% "ROI_").
% During analysis, the same folder structure is created within head_output_folder.

head_input_folder = 'D:/path/to/all/data/'
head_output_folder = 'D:/path/for/all/results/'
mkdir (head_output_folder) %create if necessary

timePoints = {"P10", "P12", "P16", "P20"}


%% settings for orientation analysis
isNematic = true;

safetyMargin = 10;  % to safely exclude influences of the boundary even at
% imprecise marking of the ROI, do not use pixels from the safetyMargin for
% orientation analysis
% This safetyMargin counts in original pixels (prior to any downsampling)
% for any resolution level

%% *FIXED* parameters of orientation analysis, assuming *ALL* images are
%% acquired with the same imaging resolution (5000-fold in the publication)
%  To compare between different imaging resolutions, vary
%  1. downsampling on the level of input grayscale image (FactordsG), 
%  2. downsampling at the level of orientations (FactordsO), and/or 
%  3. size of neighbourhood for local coherency (boxRadius_LocalAlignment)
method = 'gradient'; %for alternative method 'fourier' not used here 
                     % see https://github.com/KarlHoffmann/KarenSoans-OpticalCupAnalysis/blob/main/AnalyseOrientations/ComputeDirecty_FourierAndQ.m
boxRadius = 5; %plays the role of stencil size for method = 'gradient'
% settings for downsampling of grayscale and orientations (here implemented
% but without effect due to value 1)
FactordsG = 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
FactordsO = 1
FactordsGaO = FactordsG * FactordsO; % joint downsampling Factor Grayscale and Orientation
boxRadius_LocalAlignment = 12 % consider boxes of size side length 2*12+1=25
powerOfMagn_LocalAlignment = 1;


%% prepare to collect coherencies
filename_coherencies_collection = [head_output_folder ...
'OverviewCoherencyResults_' method '_Referenceradius' num2str(boxRadius)];

filename_coherencies_collection = [filename_coherencies_collection ...
    '_downsample' num2strForFile(FactordsG) 'x' num2strForFile(FactordsO) ...
    '_coherencyRadius' num2str(boxRadius_LocalAlignment) 'Power' num2str(powerOfMagn_LocalAlignment) '.txt'];


column_names_file = ...
        [ 'ROI #'                        ';' ...
          'boxRadius'                    ';' ...
          'boxRadius_LocalAlignment'     ';' ...
          'powerOfMagn_LocalAlignment'   ';' ...
          'global coherency'             ';' ...
          'mean local coherency'         ';' ...
          'angle of global coherency'    ';' ...
          'angle of global coherency in degree'    ';' ...
          'topfolder'                    ';' ...
          'subfolder'                    ';' ...
          'input image'                  ';' ...
          'method'                       ';' ...
          'FactordsG'                    ';' ...
          'FactordsO'                   ];

if ~exist(filename_coherencies_collection,'file')
    dlmwrite(filename_coherencies_collection, ...
        column_names_file, ...
        'delimiter','')  %header to be used
    
else
    warning(['appending coherency results without care for ordering to existing file  ' filename_coherencies_collection])
    warning('results of some images might appear multiple times')
end

coherencies_collection_arr = cell(0,14);


%% load data
for timePoint = timePoints
input_folder  = [head_input_folder char(timePoint) '/']

mkdir([head_output_folder char(timePoint)])

%% generate consistent location for .roi files that hold single ROI for one image
% Filenames of ROIs can read 'ROI_*.roi' '*_ROI.roi'
ROI_single    =  [ dir([input_folder 'ROI_*.roi']) ;  dir([input_folder '*_ROI.roi']) ] ;

for kk_ROI = 1:numel(ROI_single)
    currROI = ROI_single(kk_ROI);
    newFolderName = strrep( strrep(currROI.name, 'ROI_', 'ROIs_') , '.roi', '')
    mkdir( [currROI.folder filesep  newFolderName] )
    copyfile( [currROI.folder filesep currROI.name], [currROI.folder filesep newFolderName filesep currROI.name] )
end

% unzip ROIs in .zip files
ROI_multiple_zip =  [ dir([input_folder 'ROIs_*.zip']) ;  dir([input_folder '*_ROIs*.zip'])  ] ;


for kk_ROI = 1:numel(ROI_multiple_zip)
    ROIfolder = [ROI_multiple_zip(kk_ROI).folder filesep ROI_multiple_zip(kk_ROI).name(1:end-4) ];

    if ~exist(ROIfolder, 'dir')
        % try to unzip
        if ~exist( [ROIfolder '.zip'], 'file')
            warning('no ROI data available')
        end
        unzip( [ROIfolder '.zip'], ROIfolder)
    end
end

% get all images within the folder - adapt file format '.bmp' if needed
imagesToAnalyse  =  dir([input_folder '*.bmp']);

%% analyse all images with ROI(s) of membrane stacks
for kk_image = 1:numel(imagesToAnalyse)
    input_file2 = imagesToAnalyse(kk_image).name;
    disp(['Starting analysis of image ' input_file2])
    datetime('now')

    input_fullPath = [input_folder input_file2 ];

    ROIs_paths_struct = [ dir([imagesToAnalyse(kk_image).folder  filesep  'ROI*_' imagesToAnalyse(kk_image).name(1:end-4) '.*']) ; ...
                          dir([imagesToAnalyse(kk_image).folder  filesep  imagesToAnalyse(kk_image).name(1:end-4) '_ROI*.*'])    ] ;
        % the '.*' forces the selection to files:
        %  ROI_<filename>.roi  or  <filename>_ROI.roi  or  ROIs_<filename>.zip  or  <filename>_ROIs.zip   depending on situation (single ROI / multiple ROIs)

    for kk_prelim_path = 1: numel(ROIs_paths_struct)
        try
            if strcmp( ROIs_paths_struct(kk_prelim_path).name(end-3:end), '.zip') %can cause bad subscript errors
                % picked up a .zip file that was unpacked earlier
                % so get the individual ROIs one level deeper
                ROIs_paths_struct = dir( [ROIs_paths_struct(kk_prelim_path).folder filesep  ROIs_paths_struct(kk_prelim_path).name(1:end-4)] ) ;
                break
            end
        catch MEerror
            if ~strcmp(MEerror.identifier, 'MATLAB:badsubscript')
                MEerror.rethrow %do not silence any other error, especially not a failed assertion
            end
        end
    end

    %adapt data for which ROIs_paths_struct contains two leading entries '.' and '..'
    % from system path by removing them if present
    while numel(ROIs_paths_struct)>0 && any( strcmp(ROIs_paths_struct(1).name,  {'.', '..'} ) )
        ROIs_paths_struct(1) = [];
    end
    
    
    output_fullPath = strrep(input_fullPath, head_input_folder, head_output_folder)
    
    mkdir(fileparts(output_fullPath)) % generate folder for the output if necessary
    
    
    %% read image data
    % The name 'GrayImageStack' is used here for consistency with other analysis pipelines as 
    % in https://github.com/KarlHoffmann/KarenSoans-OpticalCupAnalysis/blob/main/AnalyseOrientations
    % Here, expect only a 2D image
    % If the data has more than one slice, consider 'for useSlice=1:NumSlices' to iterate over
    GrayImageStack_orig = imread( input_fullPath);
    [size_y_orig, size_x_orig, NumSlices] = size(GrayImageStack_orig);
    assert(NumSlices==1, ['Expecting 2D data from file ' input_fullPath ' but got ' num2str(NumSlices) ' slices.'])
    
    % normalize data
    GrayImageStack_orig = double(GrayImageStack_orig);  %convert to double if input is uint as from .bmp file
    if max(GrayImageStack_orig(:)) > 1
        GrayImageStack_orig = GrayImageStack_orig/255.0;
    end
    assert( max(GrayImageStack_orig(:)) <= 1 + 5*eps )
    
    
    %% in general: downsampling dependent on image resolution
    % for the publication, the image resolution was kept constant at 5000-fold
    
    % get image magnification from metadata
    % read .txt file with settings from TEM imaging (for the publication written by 
    % transmission electron microscopes FEI Morgagni D268 or Jeol JEM1400 Plus)
    CM_MAG = NaN; % paranoia: set to NaN to avoid re-using values from previous image
    fid = fopen( strrep(input_fullPath,'.bmp', '.txt') );
    tline = '';
    while ischar(tline) && ~startsWith(tline, '$CM_MAG ')
        tline = fgetl(fid);
    end
    CM_MAG = str2double( strrep(tline, '$CM_MAG ', '') );
    fclose(fid);
    
    % >>> Add here your downsampling factors FactordsG and FactordsO as a 
    % function of image magnification. They must be integer. <<<
    % For the publication, use fixed values of 1 but assert the image magnification
    assert(CM_MAG == 5000, ['Expecting image of magnification 5000, but got ' num2str(CM_MAG)])
    
    
    
    %% Downsampling of the grayscale image (without effect for FactordsG = 1)
    % NOTE: Identification of the foreground happens in the original image,
    % orientation analysis is done is in the downsampled grayscale data
    assert( int8(FactordsG) == FactordsG, "must be integer")
    
    GrayImageStack_blocks = 1/(FactordsG^2) ...
        * conv2in2or3D(GrayImageStack_orig, ones(FactordsG), 'valid') ;
    GrayImageStackdsG = GrayImageStack_blocks(1:FactordsG:end, 1:FactordsG:end, :);
    
    size_y_dsG = size(GrayImageStackdsG,1);
    size_x_dsG = size(GrayImageStackdsG,2);
    
    
    
    %% prepare all filenames
    %% filename to save and load back result of orientation analysis
    filename_intermediate_output = [ output_fullPath(1:end-4) '_' method '_radius' num2str(boxRadius)];

    filename_intermediate_output = [filename_intermediate_output '_downsample' num2strForFile(FactordsG) 'x' num2strForFile(FactordsO) ...
        '.mat'] ;

    %% filename to save and load back result of orientation and alignment analysis
    filename_intermediate_output1 = [filename_intermediate_output(1:end-4) ...
        '_coherencyRadius' num2str(boxRadius_LocalAlignment) 'Power' num2str(powerOfMagn_LocalAlignment) '.mat'];


    %% filename to save and load back result of orientation and alignment analysis, also per ROI
    filename_intermediate_output2 = [filename_intermediate_output1(1:end-4) ...
        'ROIwise'  '.mat'];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %% orientation analysis
    % load from pre-existing file if existent: filename_intermediate_output
    %%%%%%%%%%%%%%%%%%%%%%%
    %% try to load saved data
    canUseLoadedData = false ; %paranoia: first set to false
    for filenameToLoadFrom = {filename_intermediate_output, filename_intermediate_output1, filename_intermediate_output2}
      disp(['trying to load from file ' filenameToLoadFrom{1} ' ...'])
      if ~exist(filenameToLoadFrom{1}, 'file')
        warning('... but requested file does not exist (yet / anymore)')
      else
        try
          LoadedFromFile = load(filenameToLoadFrom{1});
          %%% FIRST CHECK that the file contains the same settings you currently have,
          %%% and the same data (as far as you already have it)
          if    (isequaln( LoadedFromFile.GrayImageStack_orig , GrayImageStack_orig ) && ...  %isequalN treats NaN as EQUAL to other NaN
                 isequal ( LoadedFromFile.input_fullPath      , input_fullPath) && ...
                 isequal ( LoadedFromFile.output_fullPath     , output_fullPath) && ...
                 isequal ( LoadedFromFile.method              , method ) && ...
                 isequal ( LoadedFromFile.boxRadius           , boxRadius ) && ...
                 isequal ( LoadedFromFile.FactordsG           , FactordsG ) && ...
                 isequal ( LoadedFromFile.FactordsO           , FactordsO )  )
            % only then allow assignment of the data you want to have
            thetadsG = LoadedFromFile.thetadsG;
            magndsG  = LoadedFromFile.magndsG;
            canUseLoadedData = true;
            disp('... was sucessful')
            break
          else
            warning('... but data in file does not fit the requested settings')
            canUseLoadedData = false;
          end
        catch %<-- probeweise auskommentieren ab hier 
          warning('... but requested file does not hold all required fields')
          canUseLoadedData = false;
        end
      end
    end

    %% analyse the orientations only if not loaded - full image in one run
    if ~canUseLoadedData
        time_needed_total=0;
        tic
        disp(['Orientation analysis by method "' method '" ...'])
        assert( isequal(method, 'gradient') && ...
                any(boxRadius == [2 3 5]), "gradient not implemented for that radius" )
        directy      = ComputeDirectyGradient(GrayImageStackdsG, boxRadius );

        thetadsG = atan2( directy(:,:,2), directy(:,:,1) );
        magndsG  = directy(:,:,3);

        time_needed=toc;
        time_needed_total = time_needed_total + time_needed;
        disp(['... done in ', num2str(time_needed), ' seconds.'])
        
        %% save newly calculated orientations to file
        % make sure not to overwrite an existing file (if file exists, we
        % should have loaded from it!)
        assert( ~exist(filename_intermediate_output,'file'), 'Orientation field was newly calculated although file with it existed!' )
        disp('Writing orientation results to file')
        
        save( filename_intermediate_output, ...
            'GrayImageStack_orig', ...    %the input data
            'input_fullPath', ...  % its source
            'output_fullPath', ...  % where data goes
            'method', 'boxRadius', ...  %the settings for orientation analysis
            'FactordsG', 'FactordsO', ... %the settings for downsampling
            'thetadsG', 'magndsG')  %the results
    end
    
    
    %% downsampling of the orientation field (without effect for FactordsO = 1)
    % NOTE that parameter  FactordsO  was set earlier
    % the input orientation field thetadsG corresponds to the first downsampling of the
    % grayscale image, and the result of the second downsampling here corresponds to
    % downsampling twice (Grayscale and Orientation - thetadsGaO)
    % as before holds
    %size_y_orig = size(GrayImageStack_orig,1);
    %size_x_orig = size(GrayImageStack_orig,2);
    %size_y_dsG = size(GrayImageStackdsG,1);
    %size_x_dsG = size(GrayImageStackdsG,2);
    
    assert( int8(FactordsO) == FactordsO, "must be integer")
    
    [thetadsGaO, magndsGaO] = downSamplDirecty_ThetaStrength(thetadsG, magndsG, FactordsO);  % skipRows=0 and skipCols=0 fits to the grayscale downsampling in this file, strength_exponent=1
    % assert theta values in the interval [0, pi] for nematic and in [0, 2*pi] for polar
    thetadsGaO =  mod(thetadsGaO, pi*(1+ not(isNematic)) );
    
    size_y_dsGaO = size(thetadsGaO,1);
    size_x_dsGaO = size(thetadsGaO,2);
    
    % do also downsample the grayscale intensities (again) as they will be used as weights
    GrayImageStack_blocksCoarse = 1/(FactordsO^2) ...
        * conv2in2or3D(GrayImageStackdsG, ones(FactordsO), 'valid') ;
    
    GrayImageStackdsGaO = GrayImageStack_blocksCoarse(1:FactordsO:end, 1:FactordsO:end, :);
    
    assert( size_y_dsGaO == size(GrayImageStackdsGaO, 1)  && ...
            size_x_dsGaO == size(GrayImageStackdsGaO, 2)  && ...
            NumSlices    == size(GrayImageStackdsGaO, 3) )  % comparison works even when GrayImageStackdsGaO is only two-dimensional




    %%%%%%%%%%%%%%%%%%%%%%%
    %% compute local alignment - potentially time consuming
    % load from pre-existing file if existent: filename_intermediate_output1
    %%%%%%%%%%%%%%%%%%%%%%%
    %% try to load saved data
    canUseLoadedData = false ; %paranoia: first set to false
    for filenameToLoadFrom = {filename_intermediate_output1, filename_intermediate_output2}
      disp(['trying to load from file ' filenameToLoadFrom{1} ' ...'])
      if ~exist(filenameToLoadFrom{1}, 'file')
        warning('... but requested file does not exist (yet / anymore)')
      else
        try
          LoadedFromFile = load(filenameToLoadFrom{1});
          %%% FIRST CHECK that the file contains the same settings you currently have,
          %%% and the same data (as far as you already have it)
          if    (isequaln( LoadedFromFile.GrayImageStack_orig , GrayImageStack_orig ) && ...  %isequalN treats NaN as EQUAL to other NaN
                 isequal ( LoadedFromFile.input_fullPath      , input_fullPath) && ...
                 isequal ( LoadedFromFile.output_fullPath     , output_fullPath) && ...
                 isequal ( LoadedFromFile.method              , method ) && ...
                 isequal ( LoadedFromFile.boxRadius           , boxRadius ) && ...
                 isequal ( LoadedFromFile.FactordsG           , FactordsG ) && ...
                 isequal ( LoadedFromFile.FactordsO           , FactordsO ) && ...
                 ... %additionally to be checked in file with local alignment data:
                 isequaln( LoadedFromFile.thetadsGaO        , thetadsGaO ) &&  ...  %isequalN treats NaN as EQUAL to other NaN
                 isequaln( LoadedFromFile.magndsGaO         , magndsGaO ) &&  ...  %isequalN treats NaN as EQUAL to other NaN
                 isequal ( LoadedFromFile.boxRadius_LocalAlignment   , boxRadius_LocalAlignment ) && ...
                 isequal ( LoadedFromFile.powerOfMagn_LocalAlignment , powerOfMagn_LocalAlignment )    )
            % only then allow assignment of the data you want to have
            localOrient     = LoadedFromFile.localOrient;
            localOrientAngle= LoadedFromFile.localOrientAngle;
            localCoherency  = LoadedFromFile.localCoherency;
            canUseLoadedData = true;
            disp('... was sucessful')
            break
          else
            warning('... but data in file does not fit the requested settings')
            canUseLoadedData = false;
          end
        catch
            warning('... but requested file does not hold all required fields')
            canUseLoadedData = false;
        end
      end
    end


    %% analyse the local alignment only if not loaded
    if ~canUseLoadedData
        % NOTE that parameters were set earlier
        disp('Calculating local alignment...')
        % try to predict runtime from analysing a small data portion
        %pick a data portion off the boundary
        yselection = max(1,min(100,size(thetadsGaO,1)-100)) : min(200,size(thetadsGaO,1)) ;
        xselection = max(1,min(100,size(thetadsGaO,2)-100)) : min(200,size(thetadsGaO,2)) ;
        tic
        ComputeLocalAlignment(thetadsGaO(yselection,xselection),magndsGaO(yselection,xselection),boxRadius_LocalAlignment, powerOfMagn_LocalAlignment);
        time_needed = toc;
        disp([ 'Expecting runtime in the range of ' ...
                num2str(time_needed / ((numel(yselection)                           ) * (numel(xselection)                           ) ) ...
                                    * ((size(thetadsGaO,1)                          ) * (size(thetadsGaO,2)                          ) ) ) ...
                ' to ' ...
                num2str(time_needed / ((numel(yselection) -2*boxRadius_LocalAlignment) * (numel(xselection) -2*boxRadius_LocalAlignment) ) ...
                                    * ((size(thetadsGaO,1)-2*boxRadius_LocalAlignment) * (size(thetadsGaO,2)-2*boxRadius_LocalAlignment) ) ) ...
                ' seconds! (based on a small test portion)'])

        %actual calculation:
        disp('Now performing local alignment analysis for full image ...')
        tic
        [localOrient, localCoherency] = ComputeLocalAlignment(thetadsGaO,magndsGaO,boxRadius_LocalAlignment, powerOfMagn_LocalAlignment);
        % give local orientation also as an angle
        localOrientAngle              = atan2(localOrient(:,:,2), localOrient(:,:,1));
        time_needed = toc;
        disp(['... done in ' num2str(time_needed) ' seconds.'])

        %% save newly calculated local alignment to file
        % make sure not to overwrite an existing file (if file exists, we
        % should have loaded from it!)
        assert( ~exist(filename_intermediate_output1,'file'), 'Local Alignment was newly calculated although file with it existed!' )
        save( filename_intermediate_output1, ...
            'GrayImageStack_orig', ...    %the input data
            'input_fullPath', ...  % its source
            'output_fullPath', ...  % where data goes
            'method', 'boxRadius', ...  %the settings for orientation analysis
            'FactordsG', 'FactordsO', ... %the settings for downsampling
            'thetadsG', 'magndsG', ...  %the results
            'thetadsGaO', 'magndsGaO', ...  %the results
            'boxRadius_LocalAlignment', 'powerOfMagn_LocalAlignment', ... %the settings for local alignment analysis
            'localOrient', 'localOrientAngle', 'localCoherency'  )    %the results

        %once the file is written, remove file which contains less information
        assert( exist(filename_intermediate_output1,'file')==2, 'Writing seems to have failed' )

        if exist(filename_intermediate_output,'file')
            delete(filename_intermediate_output)
        end
    end
    
    
    
    %% analyis per ROI
    % read ROIs if present
    
GrayImageStackROIs       = zeros(size(GrayImageStack_orig));
GrayImageStackROIs_dsGaO = zeros(size(GrayImageStackdsGaO));
for ll = 1:numel(ROIs_paths_struct)
    aROI_path = ROIs_paths_struct(ll).name;
    aROI = ReadImageJROI([ROIs_paths_struct(ll).folder filesep aROI_path]);
    try
        aROIcorners_y = aROI.vfShapes(3:3:end);
        aROIcorners_x = aROI.vfShapes(2:3:end);
        useVfShapes = true;
    catch ME
        aROIcorners_x = aROI.mnCoordinates(:,1);
        aROIcorners_y = aROI.mnCoordinates(:,2);
        useVfShapes = false;
    end
    
    %% curation of ROI data
    % in some data, aROI.vfShapes leades to inconsistent sizes of aROIcorners_x and aROIcorners_y
    % or to unrealistic large jumps in the preliminary coordinates
    % Then somewhere in between one element is missing or added, and the
    % assignment to x,y,z coordinates is mixed up
    if (length(aROIcorners_x) ~= length(aROIcorners_y)) || ...
            max( abs(aROIcorners_x(2:end) - aROIcorners_x(1:end-1)) ) > 100 || ...
            max( abs(aROIcorners_y(2:end) - aROIcorners_y(1:end-1)) ) > 100  %unrealistic large jumps
        warning('ROI might be wrong')
        if useVfShapes
            warning('trying to correct coordinates')
            CoordsCurated = aROI.vfShapes;

            %compare differences that occurs on 3-step [should be smallest]
            %to that on 2-step and 4-step
            diff2 =  CoordsCurated(3:end-1) - CoordsCurated(1:end-3)  ;
            diff3 =  CoordsCurated(4:end)   - CoordsCurated(1:end-3)  ; %the normal stepping
            diff4 = [CoordsCurated(5:end)   - CoordsCurated(1:end-4); Inf];
            diff5 = [CoordsCurated(6:end)   - CoordsCurated(1:end-5); Inf; Inf];

            IndLikelyInsertion = find( abs(diff4) < abs(diff3) );
            % if there are more than 4 mismatch indices, try to find a
            % quadruple of contiguous ones
            while numel(IndLikelyInsertion) >= 4
                warning('Possibly more than one position with extra entry in the ROI vector')
                warning('Trying to curate the first one')
                FourContiguousIndices = IndLikelyInsertion(1:end-3)+3 == IndLikelyInsertion(4:end);
                if any( FourContiguousIndices )
                    startWithinIndLikelyInsertion = find(FourContiguousIndices, 1);
                    IndLikelyInsertionFirst = IndLikelyInsertion(startWithinIndLikelyInsertion:startWithinIndLikelyInsertion+3);
                    IndToDelete = min(IndLikelyInsertionFirst)+3;
                    CoordsCurated(IndToDelete) = [];
                    warning(['removed an extra entry in the ROI coordinates vector at position ' num2str(IndToDelete) ])
                    % recompute differences after a curation step
                    diff2 =  CoordsCurated(3:end-1) - CoordsCurated(1:end-3)  ;
                    diff3 =  CoordsCurated(4:end)   - CoordsCurated(1:end-3)  ; %the normal stepping
                    diff4 = [CoordsCurated(5:end)   - CoordsCurated(1:end-4); Inf];
                    diff5 = [CoordsCurated(6:end)   - CoordsCurated(1:end-5); Inf; Inf];
                    % recompute, as indices shift behind the performed deletion
                    IndLikelyInsertion = find( abs(diff4) < abs(diff3) );
                else
                    TwoContiguousIndices = IndLikelyInsertion(1:end-1)+1 == IndLikelyInsertion(2:end);
                    if any( TwoContiguousIndices )
                        startWithinIndLikelyInsertion = find(TwoContiguousIndices, 1);
                        IndLikelyInsertion = IndLikelyInsertion(startWithinIndLikelyInsertion:startWithinIndLikelyInsertion+1);
                    end
                    break
                end
            end
            % if there are exactly 2 or exactly 4 contiguous mismatch
            % indices, they are at the deletion or insertion of a value
            if     ( numel(IndLikelyInsertion)==4 && min(IndLikelyInsertion)+3==max(IndLikelyInsertion) ) ...
                || ( numel(IndLikelyInsertion)==3 && min(IndLikelyInsertion)+2==max(IndLikelyInsertion) )
                % found a unique insertion [which is the forth of four entries in
                % IndLikelyInsertion, OR behind the third of three entries in IndLikelyInsertion]
                % and remove it
                IndToDelete = min(IndLikelyInsertion)+3;
                CoordsCurated(IndToDelete) = [];
                warning(['removed an extra entry in the ROI coordinates vector at position ' num2str(IndToDelete) ])
                % recompute differences after a curation step
                diff2 =  CoordsCurated(3:end-1) - CoordsCurated(1:end-3)  ;
                diff3 =  CoordsCurated(4:end)   - CoordsCurated(1:end-3)  ; %the normal stepping
                diff4 = [CoordsCurated(5:end)   - CoordsCurated(1:end-4); Inf];
                diff5 = [CoordsCurated(6:end)   - CoordsCurated(1:end-5); Inf; Inf];
            end

            IndLikelyDeletion1 = find( abs(diff5) < abs(diff3) );
            IndLikelyDeletion2 = find( abs(diff2) < abs(diff3) );

            if     ( numel(IndLikelyDeletion1)==3 && min(IndLikelyDeletion1)+2==max(IndLikelyDeletion1) )  ...
                && ( all(ismember(min(IndLikelyDeletion1)+[1 2], IndLikelyDeletion2) )  )
                % found a unique deletion [which is behind the third of three entries in IndLikelyDeletion1,
                % and testified by also diff2 being smaller for the second and third entry]
                % and remove it
                IndToReinsert = min(IndLikelyDeletion1)+3;
                CoordsCurated(IndToReinsert+1:end+1) = CoordsCurated(IndToReinsert:end);
                CoordsCurated(IndToReinsert) = floor(0.5*CoordsCurated(IndToReinsert-3) + 0.5*CoordsCurated(IndToReinsert+3));
                warning(['added a missing entry in the ROI coordinates vector at position ' num2str(IndToReinsert) ])
            end

            aROIcorners_y = CoordsCurated(3:3:end);
            aROIcorners_x = CoordsCurated(2:3:end);
            assert(numel(aROIcorners_y) == numel(aROIcorners_x))
        else
            aROIcorners_x = aROIcorners_x (1:min(length(aROIcorners_x), length(aROIcorners_y))) ;
            aROIcorners_y = aROIcorners_y (1:min(length(aROIcorners_x), length(aROIcorners_y))) ;
        end
    end

    % Create the binary mask image directly from corners of the ROI
    maskFromROI = roipoly(GrayImageStack_orig, aROIcorners_x, aROIcorners_y);
    %downsampling of the binary mask to the same resolution as thetadsGaO and magndsGaO
    maskFromROI_blocks = 1/( FactordsGaO^2) ...
        * conv2in2or3D(maskFromROI, ones(FactordsGaO), 'valid') ;
    maskFromROI_dsGaO  = maskFromROI_blocks(1:FactordsGaO:end, 1:FactordsGaO:end, :);

    assert( size_y_dsGaO == size(maskFromROI_dsGaO, 1)  && ...
            size_x_dsGaO == size(maskFromROI_dsGaO, 2)  && ...
            NumSlices    == size(maskFromROI_dsGaO, 3) )  % comparison works even when maskFromROI_dsGaO is only two-dimensional

    maskFromROI_dsGaO_inner = maskFromROI_dsGaO==1; %superpixels, for which all pixels are within the ROI
    maskFromROI_dsGaO_outer = maskFromROI_dsGaO >0; %superpixels, for which at least one pixel is within the ROI
    % if FactordsG==1 and FactordsO==1, then both maskFromROI_dsGaO_inner
    % and maskFromROI_dsGaO_outer are equal to maskFromROI
    % Here, use the more restrictive version maskFromROI_dsGaO_inner

    %add outline of trustworthy orientations = ROI minus boxRadius
    % further subtract boxRadius_LocalAlignment to get trustworthy coherencies
    roundUpToOdd = @(x) 1+2*ceil(0.5*x-0.5) ;
    maskFromROI_orientReliab = imerode( maskFromROI_dsGaO_inner,  strel('rectangle',[roundUpToOdd( 2*(boxRadius/FactordsO + safetyMargin/FactordsGaO) ), roundUpToOdd( 2*(boxRadius/FactordsO + safetyMargin/FactordsGaO) ) ]) ) ;
    maskFromROI_coherReliab  = imerode( maskFromROI_orientReliab, strel('rectangle',[2*boxRadius_LocalAlignment+1, 2*boxRadius_LocalAlignment+1]) ) ;

    thetadsGaO_masked_inner             = reshape( thetadsGaO      (maskFromROI_dsGaO_inner),  [], 1); %select and flatten
    magndsGaO_masked_inner              = reshape( magndsGaO       (maskFromROI_dsGaO_inner),  [], 1);
    thetadsGaO_masked_orientReliab      = reshape( thetadsGaO      (maskFromROI_orientReliab), [], 1);
    magndsGaO_masked_orientReliab       = reshape( magndsGaO       (maskFromROI_orientReliab), [], 1);
    localCoherency_masked_inner         = reshape( localCoherency  (maskFromROI_dsGaO_inner),  [], 1);
    localCoherency_masked_coherReliab   = reshape( localCoherency  (maskFromROI_coherReliab),  [], 1);
    localOrientAngle_masked_inner       = reshape( localOrientAngle(maskFromROI_dsGaO_inner),  [], 1);
    localOrientAngle_masked_coherReliab = reshape( localOrientAngle(maskFromROI_coherReliab),  [], 1);

    [overallOrient,                overallCoherency]                = ComputeQTensor(thetadsGaO_masked_inner,        magndsGaO_masked_inner,        powerOfMagn_LocalAlignment);
    [overallOrient_orientReliab,   overallCoherency_orientReliab]   = ComputeQTensor(thetadsGaO_masked_orientReliab, magndsGaO_masked_orientReliab, powerOfMagn_LocalAlignment);

    meanLocalCoherency             =  nanmean(localCoherency_masked_inner(:));
    meanLocalCoherency_coherReliab =  nanmean(localCoherency_masked_coherReliab(:));

    % calculation of angle of global coherency
    overallOrientAngle              = atan2(overallOrient(2), overallOrient(1));
    overallOrientAngle_orientReliab = atan2(overallOrient_orientReliab(2), overallOrient_orientReliab(1));


    % add to struct of results
    ROIs_paths_struct(ll).corners_x                                = aROIcorners_x;
    ROIs_paths_struct(ll).corners_y                                = aROIcorners_y;
    ROIs_paths_struct(ll).overallOrient                            = overallOrient;
    ROIs_paths_struct(ll).overallOrient_orientReliab               = overallOrient_orientReliab;
    ROIs_paths_struct(ll).overallCoherency                         = overallCoherency;
    ROIs_paths_struct(ll).overallCoherency_orientReliab            = overallCoherency_orientReliab;
    ROIs_paths_struct(ll).overallOrientAngle                       = overallOrientAngle;
    ROIs_paths_struct(ll).overallOrientAngle_orientReliab          = overallOrientAngle_orientReliab;
    ROIs_paths_struct(ll).overallOrientAngle_inDegree              = mod(rad2deg(overallOrientAngle), 180);
    ROIs_paths_struct(ll).overallOrientAngle_orientReliab_inDegree = mod(rad2deg(overallOrientAngle_orientReliab), 180);
    ROIs_paths_struct(ll).meanLocalCoherency                       = meanLocalCoherency;
    ROIs_paths_struct(ll).meanLocalCoherency_coherReliab           = meanLocalCoherency_coherReliab;
    ROIs_paths_struct(ll).thetadsGaO_masked_inner                  = thetadsGaO_masked_inner;
    ROIs_paths_struct(ll).magndsGaO_masked_inner                   = magndsGaO_masked_inner;
    ROIs_paths_struct(ll).thetadsGaO_masked_orientReliab           = thetadsGaO_masked_orientReliab;
    ROIs_paths_struct(ll).magndsGaO_masked_orientReliab            = magndsGaO_masked_orientReliab;
    ROIs_paths_struct(ll).localCoherency_masked_inner              = localCoherency_masked_inner;
    ROIs_paths_struct(ll).localCoherency_masked_coherReliab        = localCoherency_masked_coherReliab;
    ROIs_paths_struct(ll).localOrientAngle_masked_inner            = localOrientAngle_masked_inner;
    ROIs_paths_struct(ll).localOrientAngle_masked_coherReliab      = localOrientAngle_masked_coherReliab;
    % the pairs (localCoherency_masked_inner,       localOrientAngle_masked_inner )
    %       and (localCoherency_masked_coherReliab, localOrientAngle_masked_coherReliab )
    % act like  (magndsGaO_masked_inner, thetadsGaO_masked_inner)
    % as descriptors of some vectors (but all these without information on the relative locations of the vectors)
    % ROIs_paths_struct(ll).overallCoherency, ROIs_paths_struct(ll).meanLocalCoherency
    
    %% build up image with only foreground regions
    GrayImageStackROIs       = GrayImageStackROIs       + GrayImageStack_orig .* maskFromROI;
    GrayImageStackROIs_dsGaO = GrayImageStackROIs_dsGaO + GrayImageStackdsGaO .* maskFromROI_dsGaO_inner;

    %% plot into overview figure with all ROIs of the image
    if ~exist('FigOverview','var') || ~ishandle(FigOverview)
        %create only if not existent yet
        FigOverview = figure(1);

        subplot(1,2,1)
        imshow(GrayImageStack_orig, 'InitialMagnification', 'fit')
        hold on
    else
        % make existing figure current
        figure(FigOverview)
    end
    subplot(1,2,1)
    hold on
    plot( [aROIcorners_x ; aROIcorners_x(1)] ,...  %close the loop
          [aROIcorners_y ; aROIcorners_y(1)], 'yellow', 'LineWidth', 3)
    subplot(1,2,2)
    imshow(GrayImageStackROIs, 'InitialMagnification', 'fit')
    % alternatively to highlight ROIs:  imshow(maskFromROI, 'InitialMagnification', 'fit')
    drawnow()
    pause(0.1)
    
    % save the figure when all ROIs are enclosed
    if ll == numel(ROIs_paths_struct)
        print(FigOverview, [output_fullPath(1:end-4) '_overviewROIs'],'-dpdf', '-r0', '-fillpage')
        print(FigOverview, [output_fullPath(1:end-4) '_overviewROIs'],'-dpng', '-r0')
        pause(0.1)
        close(FigOverview)
        clear FigOverview  % also delete variable for clean build-up for the next image
    end
    
    %% plot distributions of orientations and locally averaged orientations within the ROI
    HistBinLimits = [0, 0.2];
    MyNumBins = 100;
    colOrient      = [0   0 1  ];
    colOrient_magn = [0   0 0.5];
    colReliab      = [1   0 0  ];
    colReliab_magn = [0.5 0 0  ];
    colGlobalCoh   = [1   1 0  ];

    NbinsPolar = 128;
    EdgesPolar = linspace(0, pi, NbinsPolar+1);
    % this ensures to have all data fall into the bins (NaN bins are only due to theta value NaN)
    assert ( EdgesPolar(1) == 0  &  EdgesPolar(end) == pi , 'Incorrect bins')
    EdgesNematic = [EdgesPolar(1:end-1)  pi+EdgesPolar];


    % figure layout
    FigDistr = figure('PaperOrientation', 'landscape', 'Position', [0, 0, 1400, 700]) ;
    % explicitly create the axes at specified locations [left bottom width heigth]
    FigDistrA      =      axes('Units', 'pixels', 'Position', [   0 350  650 300] + [30 20 -20 -20]) ;
    FigDistrB      = polaraxes('Units', 'pixels', 'Position', [   0   0  650 300] + [30 20 -20 -20]) ;
    FigDistrC      =      axes('Units', 'pixels', 'Position', [ 700 350  650 300] + [30 20 -20 -20]) ;
    FigDistrD      = polaraxes('Units', 'pixels', 'Position', [ 700   0  650 300] + [30 20 -20 -20]) ;
    FigDistrtitle  =      axes('Units', 'pixels', 'OuterPosition', [   0 650 1400 50]) ;

    %overall title
    axes(FigDistrtitle)
    titletext = text(0.5,1, {'empirical probability density function of orientation magnitudes', ...
               ['downsampling grayscale = ' num2str(FactordsG) ', downsampling orientations = ' num2str(FactordsO) ...
               ', local averaging box = ' num2str(2*boxRadius_LocalAlignment+1) 'x' num2str(2*boxRadius_LocalAlignment+1) ] }, ...
               'FontSize', 15, 'FontWeight', 'bold', ...
               'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    %clear coordinate axis
    FigDistrtitle.XAxis.Visible = 'off';
    FigDistrtitle.YAxis.Visible = 'off';

    %% magnitudes of orientation vectors
    axes(FigDistrA)
    title('Orientation magnitudes')
    hist_magndsGaO = histogram(FigDistrA, magndsGaO_masked_inner, ...
        'BinLimits', HistBinLimits, ...
        'NumBins', MyNumBins, ...
        'FaceColor', colOrient, ...
        'FaceAlpha', 0.5, ...
        'Normalization','pdf', ...
        'DisplayName', ['Orientation magnitudes in ROI ' num2str(ll)] );
    hold on
    hist_magndsGaO_orientReliab = histogram(FigDistrA, magndsGaO_masked_orientReliab, ...
        'BinLimits', HistBinLimits, ...
        'NumBins', MyNumBins, ...
        'FaceColor', colReliab, ...
        'FaceAlpha', 0.5, ...
        'Normalization','pdf', ...
        'DisplayName', ['Orientation magnitudes in ROI ' num2str(ll) ' where box fits'] );

    ylim_curr = ylim();
    ylim_curr(2) = max( ylim_curr(2), ... % existing value, or round to next number (1...9)*10^n
        10^floor(log10(max(hist_magndsGaO.Values))) * ...
        ceil( max(hist_magndsGaO.Values) / (10^floor(log10(max(hist_magndsGaO.Values)))) ) );
    ylim(ylim_curr)

    %show mean value
    line_magndsGaO_mean = plot( [1 1]*mean(magndsGaO_masked_inner), [0 ylim_curr(2)], ...
        'Color', colOrient, ...
        'LineWidth', 3, ...
        'DisplayName', ['Mean orientation magnitudes in ROI ' num2str(ll)] );
    line_magndsGaO_mean_orientReliab = plot( [1 1]*mean(magndsGaO_masked_orientReliab), [0 ylim_curr(2)], ...
        'Color', colReliab, ...
        'LineWidth', 3, ...
        'DisplayName', ['Mean orientation magnitudes in ROI ' num2str(ll) ' where box fits'] );
    line_coherency      = plot( [1 1]*overallCoherency, [0 ylim_curr(2)], ...
        'Color', colGlobalCoh, ...
        'LineWidth', 3, ...
        'DisplayName', ['Global coherency in ROI ' num2str(ll) ' = ' num2str(overallCoherency) ] );

    xlim([0, 0.2])
    grid on

    lgd = legend([hist_magndsGaO,     hist_magndsGaO_orientReliab, ...
                 line_magndsGaO_mean, line_magndsGaO_mean_orientReliab, ...
                 line_coherency], 'Location', 'northoutside', ...
                 'Box', 'on');
    %%% NOTE: for the publication, plotting was restricted to  hist_magndsGaO_orientReliab, 
    %%% line_magndsGaO_mean_orientReliab, and line_coherency (in orange)
    
        
    %% orientation angles of orientation vectors (optional with magnitude as weight)
    axes(FigDistrB)
    %%% calculate frequencies of orientations
    % bin the data yourself to allow for weights being considered
    % bin into interval [0,pi], then duplicate data from theta to pi+theta to represent nematics
    theta_for_hists              = mod(thetadsGaO_masked_inner, pi);
    assert ( isempty(theta_for_hists) || ...
             min(theta_for_hists(:))>= 0  &  max(theta_for_hists(:))<=pi, 'theta falls  out of range [0, pi]' )
    theta_for_hists_orientReliab = mod(thetadsGaO_masked_orientReliab, pi);
    assert ( isempty(theta_for_hists_orientReliab) || ...
             min(theta_for_hists_orientReliab(:))>= 0  &  max(theta_for_hists_orientReliab(:))<=pi, 'theta falls  out of range [0, pi]' )

    %the numbers per bin - for data of "masked_inner"
    whichBin = discretize(theta_for_hists, EdgesPolar);
    Freq                    = zeros([1 numel(EdgesPolar)-1]);
    Freq_magn               = zeros([1 numel(EdgesPolar)-1]);
    for kk=1:numel(EdgesPolar)-1
        Freq(kk)                      =    sum(   sum( (whichBin==kk) ));
        %Freq_magn(kk)                 = nansum(nansum( (whichBin==kk) .* magndsGaO_masked_inner ));
        Freq_magn(kk)                 =    sum(   sum(  magndsGaO_masked_inner(whichBin==kk) ));
    end
    
    %the numbers per bin - for data of "orientReliab"
    whichBin = discretize(theta_for_hists_orientReliab, EdgesPolar);
    Freq_orientReliab       = zeros([1 numel(EdgesPolar)-1]);
    Freq_orientReliab_magn  = zeros([1 numel(EdgesPolar)-1]);
    for kk=1:numel(EdgesPolar)-1
        Freq_orientReliab(kk)         =    sum(   sum( (whichBin==kk) ));
        Freq_orientReliab_magn(kk)    = nansum(nansum( (whichBin==kk) .* magndsGaO_masked_orientReliab ));
    end

    pHist                  = polarhistogram('BinEdges',EdgesNematic,'BinCounts',[Freq                    Freq            ], ...
        'DisplayStyle','stairs', 'EdgeColor', colOrient ...
        );  ... %not with style 'stairs':   'FaceAlpha', PolarAlpha, 'FaceColor', colOrient);
    hold on
    pHist_magn             = polarhistogram('BinEdges',EdgesNematic,'BinCounts',[Freq_magn               Freq_magn       ], ...
        'DisplayStyle','stairs', 'EdgeColor', colOrient_magn ...
        );  ... %not with style 'stairs':   'FaceColor', colOrient_magn);
    pHist_orientReliab     = polarhistogram('BinEdges',EdgesNematic,'BinCounts',[Freq_orientReliab       Freq_orientReliab     ], ...
        'DisplayStyle','stairs', 'EdgeColor', colReliab ...
        );  ... %not with style 'stairs':   'FaceColor', colReliab);
    pHist_orientReliab_magn= polarhistogram('BinEdges',EdgesNematic,'BinCounts',[Freq_orientReliab_magn  Freq_orientReliab_magn], ...
        'DisplayStyle','stairs', 'EdgeColor', colReliab_magn ...
        );  ... %not with style 'stairs':   'FaceColor', colReliab_magn);

    % add global average orientation
    pAverageOrient = polarplot( [0 pi]+atan2(overallOrient(2),overallOrient(1)), [0.8, 0.8]* max(Freq), ...
        'Color', colGlobalCoh, 'LineWidth', 3 );

    lgd = legend([pHist, pHist_magn, pHist_orientReliab, pHist_orientReliab_magn, ...
                  pAverageOrient], ...
                 {'count', 'magnitude-weighted', 'count where box fits', 'magnitude-weighted where box fits', ...
                  'global average'}, ...
        'Location', 'eastoutside');
    lgd.Title.String = ['orientations in ROI ' num2str(ll)];
    
    %%% NOTE: for the publication, plotting was restricted to  pHist_orientReliab, 
    %%% and pAverageOrient (in orange)
    
    
    %% magnitudes of locally averaged orientation vectors
    axes(FigDistrC)
    title('Local Coherency = magnitude of locally averaged orientations')
    hist_localCoh = histogram(FigDistrC, localCoherency_masked_inner, ...
        'BinLimits', [0, 1], ...
        'NumBins', 5*MyNumBins, ...
        'FaceColor', colOrient, ...
        'FaceAlpha', 0.5, ...
        'EdgeColor', colOrient, ...
        'Normalization','pdf', ...
        'DisplayName', ['Local coherencies in ROI ' num2str(ll)] );
    hold on
    hist_localCoh_orientReliab = histogram(FigDistrC, localCoherency_masked_coherReliab, ...
        'BinLimits', [0, 1], ...
        'NumBins', 5*MyNumBins, ...
        'FaceColor', colReliab, ...
        'FaceAlpha', 0.5, ...
        'EdgeColor', colReliab, ...
        'Normalization','pdf', ...
        'DisplayName', ['Local coherencies in ROI ' num2str(ll) ' where box fits'] );
    
    ylim_curr = ylim();
    ylim_curr(2) = max( ylim_curr(2), ... % existing value, or round to next number (1...9)*10^n
        10^floor(log10(max(hist_localCoh.Values))) * ...
        ceil( max(hist_localCoh.Values) / (10^floor(log10(max(hist_localCoh.Values)))) ) );
    ylim(ylim_curr)
    
    %show mean value
    line_localCoh_mean = plot( [1 1]*mean(localCoherency_masked_inner), [0 ylim_curr(2)], ...
        'Color', colOrient, ...
        'LineWidth', 3, ...
        'DisplayName', ['Mean local coherency in ROI ' num2str(ll)] );
    line_localCoh_mean_orientReliab = plot( [1 1]*mean(localCoherency_masked_coherReliab), [0 ylim_curr(2)], ...
        'Color', colReliab, ...
        'LineWidth', 3, ...
        'DisplayName', ['Mean local coherency in ROI ' num2str(ll) ' where box fits'] );
    line_coherency      = plot( [1 1]*overallCoherency, [0 ylim_curr(2)], ...
        'Color', colGlobalCoh, ...
        'LineWidth', 3, ...
        'DisplayName', ['Global coherency in ROI ' num2str(ll) ' = ' num2str(overallCoherency) ] );
    
    xlim([0, 1])
    grid on
    
    lgd = legend([hist_localCoh,     hist_localCoh_orientReliab, ...
                 line_localCoh_mean, line_localCoh_mean_orientReliab, ...
                 line_coherency], 'Location', 'northoutside', ...
                 'Box', 'on');
    
    
    %% orientation angles of locally averaged orientation vectors (optional with magnitude as weight)
    axes(FigDistrD)
    %%% calculate frequencies of orientations
    % bin the data yourself to allow for weights being considered
    % bin into interval [0,pi], then duplicate data from theta to pi+theta to represent nematics
    thetaLocalOrient_for_hists              = mod(localOrientAngle_masked_inner, pi);
    assert ( isempty(thetaLocalOrient_for_hists) || ...
             min(thetaLocalOrient_for_hists(:))>= 0  &  max(thetaLocalOrient_for_hists(:))<=pi, 'theta falls  out of range [0, pi]' )
    thetaLocalOrient_for_hists_coherReliab = mod(localOrientAngle_masked_coherReliab, pi);
    assert ( isempty(localOrientAngle_masked_coherReliab) || ...
             (min(thetaLocalOrient_for_hists_coherReliab(:))>= 0  &  max(thetaLocalOrient_for_hists_coherReliab(:))<=pi), ...
             'theta falls  out of range [0, pi]' )
    
    %the numbers per bin - for data of "masked_inner"
    whichBin = discretize(thetaLocalOrient_for_hists, EdgesPolar);
    FreqLocalCoh                    = zeros([1 numel(EdgesPolar)-1]);
    FreqLocalCoh_magn               = zeros([1 numel(EdgesPolar)-1]);
    for kk=1:numel(EdgesPolar)-1
        FreqLocalCoh(kk)                      =    sum(   sum( (whichBin==kk) ));
       %FreqLocalCoh_magn(kk)                 = nansum(nansum( (whichBin==kk) .* localCoherency_masked_inner ));
        FreqLocalCoh_magn(kk)                 =    sum(   sum(  localCoherency_masked_inner(whichBin==kk) ));
    end
    
    %the numbers per bin - for data of "orientReliab"
    whichBin = discretize(thetaLocalOrient_for_hists_coherReliab, EdgesPolar);
    FreqLocalCoh_orientReliab       = zeros([1 numel(EdgesPolar)-1]);
    FreqLocalCoh_orientReliab_magn  = zeros([1 numel(EdgesPolar)-1]);
    for kk=1:numel(EdgesPolar)-1
        FreqLocalCoh_orientReliab(kk)         =    sum(   sum( (whichBin==kk) ));
        FreqLocalCoh_orientReliab_magn(kk)    =    sum(   sum( (whichBin==kk) .* localCoherency_masked_coherReliab ));
    end
    
    pHistlocalCoh                  = polarhistogram('BinEdges',EdgesNematic,'BinCounts',[FreqLocalCoh                    FreqLocalCoh            ], ...
        'DisplayStyle','stairs', 'EdgeColor', colOrient ...
        );  ... %not with style 'stairs':   'FaceAlpha', PolarAlpha, 'FaceColor', colOrient);
    hold on
    pHistlocalCoh_magn             = polarhistogram('BinEdges',EdgesNematic,'BinCounts',[FreqLocalCoh_magn               FreqLocalCoh_magn       ], ...
        'DisplayStyle','stairs', 'EdgeColor', colOrient_magn ...
        );  ... %not with style 'stairs':   'FaceColor', colOrient_magn);
    pHistlocalCoh_orientReliab     = polarhistogram('BinEdges',EdgesNematic,'BinCounts',[FreqLocalCoh_orientReliab       FreqLocalCoh_orientReliab     ], ...
        'DisplayStyle','stairs', 'EdgeColor', colReliab ...
        );  ... %not with style 'stairs':   'FaceColor', colReliab);
    pHistlocalCoh_orientReliab_magn= polarhistogram('BinEdges',EdgesNematic,'BinCounts',[FreqLocalCoh_orientReliab_magn  FreqLocalCoh_orientReliab_magn], ...
        'DisplayStyle','stairs', 'EdgeColor', colReliab_magn ...
        );  ... %not with style 'stairs':   'FaceColor', colReliab_magn);
    
    % add global average orientation
    pAverageOrient = polarplot( [0 pi]+atan2(overallOrient(2),overallOrient(1)), [0.8, 0.8]* max(FreqLocalCoh), ...
        'Color', colGlobalCoh, 'LineWidth', 3 );
    
    lgd = legend([pHistlocalCoh, pHistlocalCoh_magn, pHistlocalCoh_orientReliab, pHistlocalCoh_orientReliab_magn, ...
                  pAverageOrient], ...
                 {'count', 'magnitude-weighted', 'count where box fits', 'magnitude-weighted where box fits', ...
                  'global average'}, ...
        'Location', 'eastoutside');
    lgd.Title.String = ['Locally averaged orientations in ROI ' num2str(ll)];
    
    
    %% save plot (into one subfolder per ROI)
    subfolderThisROI = [ output_fullPath(1:end-4) '_ROI' num2str(ll)];
    mkdir( subfolderThisROI )
    
    [~, filename_processingInfo, ~] = fileparts(filename_intermediate_output1);
    
    filename_print_singleROIdistr = [subfolderThisROI  filesep  filename_processingInfo ...
        '_ROI' num2str(ll) '_' ... %%% ONLY THE SUBFOLDER AND THE ADDED ROI NUMBER differ from filename of full image plot
        'Distributions' ] ;
    
    
    print(FigDistr, filename_print_singleROIdistr,'-dpdf', '-r0', '-fillpage')
    print(FigDistr, filename_print_singleROIdistr,'-dpng', '-r0')
    
    pause(0.1)
    close(FigDistr)
end




%% save result arrays of local alignment (together with their settings) to file
%% NOTE: the computation of the saved data can take quite a long time - so
%  CHECK TWICE whether you want to overwrite it
save( filename_intermediate_output2, ...
        'GrayImageStack_orig', ...    %the input data
        'input_fullPath', ...  % its source
        'output_fullPath', ...  % where data goes
        'method', 'boxRadius', ...  %the settings for orientation analysis
        'FactordsG', 'FactordsO', ... %the settings for downsampling
        'thetadsG', 'magndsG', ...  %the results
        'thetadsGaO', 'magndsGaO', ...  %the results
        'boxRadius_LocalAlignment', 'powerOfMagn_LocalAlignment', ... %the settings for local alignment analysis
        'localOrient', 'localOrientAngle', 'localCoherency', ...   %the results
        'ROIs_paths_struct'  )   %the information on ROIs

%once the file is written, remove file which contains less information
assert( exist(filename_intermediate_output2,'file')==2, 'Writing seems to have failed' )

if exist(filename_intermediate_output1,'file')
    delete(filename_intermediate_output1)
end



%% collect and output and the obtained coherency values (to be saved for all images together)
for ll = 1:numel(ROIs_paths_struct)
        % for output to command window
        disp(input_fullPath)
        disp(         {['ROI #' num2str(ll) ], ...
                       ['boxRadius=' num2str(boxRadius) ', boxRadiusAlignment=' num2str(boxRadius_LocalAlignment)]} )

        disp( {'global coherency', [num2str(ROIs_paths_struct(ll).overallCoherency)    ]} )
        disp( {'mean local coherency', [num2str(ROIs_paths_struct(ll).meanLocalCoherency)    ]} )
        disp( {'angle of global coherency in degree', [num2str(ROIs_paths_struct(ll).overallOrientAngle_inDegree)    ]} )

        %% collect coherencies in file
            input_PathPieces = strsplit(input_fullPath, '/');

            % collect coherencies in a cell array
            coherencies_collection_arr(end+1,:) = {ll, boxRadius, boxRadius_LocalAlignment, powerOfMagn_LocalAlignment, ROIs_paths_struct(ll).overallCoherency, ROIs_paths_struct(ll).meanLocalCoherency, ROIs_paths_struct(ll).overallOrientAngle, ROIs_paths_struct(ll).overallOrientAngle_inDegree, ...
                    input_PathPieces{end-2}, input_PathPieces{end-1}, input_PathPieces{end},  ...
                    method, FactordsG, FactordsO  };

            toBeWritten = cellfun( (@(strORnum) [num2str(strORnum) ';']) , coherencies_collection_arr(end,:), 'UniformOutput', false);
 
            % header was written outside the loops
            dlmwrite( filename_coherencies_collection , ...
                ' ', ... % force (additional) newline in the output file
                '-append', 'delimiter','', 'newline','pc')
            dlmwrite( filename_coherencies_collection , ...
                cell2mat(toBeWritten), ...
                '-append', 'delimiter','', 'newline','pc')
end



%% plot the full input image together with analysis results, 
% using different selections of which results to show

%% settings for plotting
% colors
col_director    = [0.4 0.4 1.0]; %NOT 'blue', as this is too dark
col_localOrient = [1.0 0.4 0.4]; %NOT 'red'

useMagnForPlotting = true;
centerQuiverOnPixels = true;

% ensure roughly the same spacing in physical units across different resolutions - p.e. use spacing = boxRadius;
spacingInOrigPixels = 16;
spacing = max( 1, round(spacingInOrigPixels / FactordsGaO) );


plotSettings = struct( ...
    ... % each "column" below corresponds to one plot setting
    ... % expand at your dicretion 
    'showLocalCoherencyInsteadOfOrig', [false, true,  false], ...
    'showOrient',                      [true,  false, false], ...
    'show_localOrient',                [false, true,  false], ...
    'write_ROICoherency',              [false, false, true ], ...  %write coherency values on the plot
    'show_ROIwiseOrient',              [false, true,  false] ...
    )


for iteratorOverPlotSettings = 1:numel(plotSettings.showOrient)
    showLocalCoherencyInsteadOfOrig = plotSettings.showLocalCoherencyInsteadOfOrig(iteratorOverPlotSettings);
    showOrient                      = plotSettings.showOrient                     (iteratorOverPlotSettings);
    show_localOrient                = plotSettings.show_localOrient               (iteratorOverPlotSettings);
    write_ROICoherency              = plotSettings.write_ROICoherency             (iteratorOverPlotSettings);
    show_ROIwiseOrient              = plotSettings.show_ROIwiseOrient             (iteratorOverPlotSettings);

    %% prepare orientation field and averaged local orientation field
    u = magndsGaO.^useMagnForPlotting .*  cos( thetadsGaO ) ; %    .* GrayImageStack;
    v = magndsGaO.^useMagnForPlotting .* -sin( thetadsGaO ) ; %    .* GrayImageStack; %negative sign to transform from theta in normal coord layout to coord layout used by MATLAB (y downwards in plotting that is started by imshow() )
    u_localOrient = localOrient(:,:,1) .* localCoherency;
    v_localOrient = localOrient(:,:,2) .* localCoherency;
    [x12, y12] = meshgrid(0.5+(0.5+ 0:size(GrayImageStackdsGaO,2) )*FactordsGaO, ...
                         0.5+(0.5+ 0:size(GrayImageStackdsGaO,1) )*FactordsGaO );
    % omit last entry in x12/y12 if their size exeeds size of u/v (is the case when downsampling had exactly fitting blocks)
    x12 = x12(1:size(u,1), 1:size(u,2));
    y12 = y12(1:size(u,1), 1:size(u,2));
    
    %% plotting
    Fig12 = figure('PaperOrientation', 'landscape', 'Position', [0, 0, 1200, 700]);
    [axesHandles12, positions12] = tight_subplot(1,1,[.08 .05],[.01 .05],[.05 .02]) ;   %gap, margin_height, margin_width)
    fig12A = axesHandles12(1);
    
    axes(fig12A)   %make the existing axes current
    
    %%% original grayscale OR local coherency
    if showLocalCoherencyInsteadOfOrig
        OrigData12 = imshow(localCoherency, 'XData', [0.5, 0.5+size(GrayImageStack_orig,2)], ...
                                            'YData', [0.5, 0.5+size(GrayImageStack_orig,1)]    );
    else
        OrigData12 = imshow(GrayImageStack_orig);
    end
    hold on
    
    % quiver for orientations
    if showOrient
      quiv_onOrigData1 = quiver(x12(1:spacing:end,1:spacing:end), y12(1:spacing:end,1:spacing:end), u(1:spacing:end,1:spacing:end), v(1:spacing:end,1:spacing:end));
      if centerQuiverOnPixels
        %duplicate quiver (with opposite sign) to have coherent appearance
        quiv_onOrigData2 = quiver(x12(1:spacing:end,1:spacing:end), y12(1:spacing:end,1:spacing:end), -u(1:spacing:end,1:spacing:end), -v(1:spacing:end,1:spacing:end));
      end

      % modify quiver plot(s)
      quiv_onOrigData1.ShowArrowHead = 'off';
      quiv_onOrigData1.LineWidth  = 1.5;  %(default: 0.5)
      quiv_onOrigData1.AutoScaleFactor = 2.0 ;
      quiv_onOrigData1.Color = col_director;
      if centerQuiverOnPixels
        %give the duplicate quiver the same properties, especially the same color
        quiv_onOrigData2.ShowArrowHead   = quiv_onOrigData1.ShowArrowHead   ;
        quiv_onOrigData2.LineWidth       = quiv_onOrigData1.LineWidth       ;
        quiv_onOrigData2.AutoScaleFactor = quiv_onOrigData1.AutoScaleFactor ;
        quiv_onOrigData2.Color           = quiv_onOrigData1.Color           ;
      end
    end
    
    
    if show_localOrient
      quiv_localOrient1 = quiver(x12(1:spacing:end,1:spacing:end), y12(1:spacing:end,1:spacing:end), u_localOrient(1:spacing:end,1:spacing:end), -v_localOrient(1:spacing:end,1:spacing:end));
      %Note negative sign on y-component to transform from theta in normal coord layout to coord layout used by MATLAB (y downwards in plotting that is started by imshow() )
      if centerQuiverOnPixels
        %duplicate quiver (with opposite sign) to have coherent appearance
        quiv_localOrient2 = quiver(x12(1:spacing:end,1:spacing:end), y12(1:spacing:end,1:spacing:end), -u_localOrient(1:spacing:end,1:spacing:end), v_localOrient(1:spacing:end,1:spacing:end));
      end
      
      % modify quiver plot(s)
      quiv_localOrient1.ShowArrowHead = 'off';
      quiv_localOrient1.LineWidth  = 0.5;  %(default: 0.5)
      quiv_localOrient1.AutoScaleFactor = 1.0 ;
      quiv_localOrient1.Color = col_localOrient;
      if centerQuiverOnPixels
        %give the duplicate quiver the same properties, especially the same color
        quiv_localOrient2.ShowArrowHead   = quiv_localOrient1.ShowArrowHead   ;
        quiv_localOrient2.LineWidth       = quiv_localOrient1.LineWidth       ;
        quiv_localOrient2.AutoScaleFactor = quiv_localOrient1.AutoScaleFactor ;
        quiv_localOrient2.Color           = quiv_localOrient1.Color           ;
      end
    end
    
    
    
    % add outlines of ROIs
    for ll = 1:numel(ROIs_paths_struct)
        % original ROIs
        plot( [ROIs_paths_struct(ll).corners_x ; ROIs_paths_struct(ll).corners_x(1)] ,...  %close the loop
              [ROIs_paths_struct(ll).corners_y ; ROIs_paths_struct(ll).corners_y(1)], 'yellow', 'LineWidth', 3)
        
        %add outline of trustworthy orientations = ROI minus boxRadius (here calculated in original pixel numbers, therefore FactordsG added)
        % further subtract boxRadius_LocalAlignment to get trustworthy coherencies (here calculated in original pixel numbers, therefore FactordsGaO added)
        maskFromROI = roipoly(GrayImageStack_orig, ROIs_paths_struct(ll).corners_x, ROIs_paths_struct(ll).corners_y);
        maskFromROI_orientReliab = imerode( maskFromROI,              strel('rectangle',FactordsG*[2*boxRadius+1,2*boxRadius+1]) ) ;
        maskFromROI_coherReliab  = imerode( maskFromROI_orientReliab, strel('rectangle',FactordsGaO*[2*boxRadius_LocalAlignment+1,2*boxRadius_LocalAlignment+1]) ) ;
        
        [boundaries, labels] = bwboundaries(maskFromROI_orientReliab);
        for bb = 1:length(boundaries)
            boundary = boundaries{bb};
            plot(boundary(:,2), boundary(:,1), 'yellow', 'LineStyle', ':', 'LineWidth', 2)
        end
        
        [boundaries, labels] = bwboundaries(maskFromROI_coherReliab);
        for bb = 1:length(boundaries)
            boundary = boundaries{bb};
            plot(boundary(:,2), boundary(:,1), 'yellow', 'LineStyle', ':', 'LineWidth', 2)
        end
        
        
        if write_ROICoherency % writing into the plot
            text( -100+mean(ROIs_paths_struct(ll).corners_x), -100+mean(ROIs_paths_struct(ll).corners_y), ...
                {'global coherency', [num2str(ROIs_paths_struct(ll).overallCoherency)    ]});%    ' / ' num2str(ROIs_paths_struct(ll).overallCoherency_orientReliab) ]})
            
            text( -100+mean(ROIs_paths_struct(ll).corners_x), +100+mean(ROIs_paths_struct(ll).corners_y), ...
                {'mean local coherency', [num2str(ROIs_paths_struct(ll).meanLocalCoherency)    ]});%   ' / ' num2str(ROIs_paths_struct(ll).meanLocalCoherency_coherReliab) ]})
        end
        
        
        if show_ROIwiseOrient
            quiver( mean(ROIs_paths_struct(ll).corners_x)*[1 1],                      mean(ROIs_paths_struct(ll).corners_y)*[1 1], ...
                        (ROIs_paths_struct(ll).overallOrient_orientReliab(1))*[1 -1],     (ROIs_paths_struct(ll).overallOrient_orientReliab(2))*[-1 +1], ...
                    'red', 'ShowArrowHead', 'off',...
                    'LineWidth', 5, ...  %(default: 0.5)
                    'AutoScaleFactor', 50)
        end
        
        % finalise image - only inside the loop as output focussed on each ROI is wanted
        axis equal
        title({[repmat('Coherency field', showLocalCoherencyInsteadOfOrig) ...
            repmat('orig image', 1-showLocalCoherencyInsteadOfOrig) ...
            ' with ROI #' num2str(ll) ' outlined, ' ...  %%% ONLY THIS LINE DIFFERS from title of full image plot
            repmat(['orientations (blue, shown every ' num2str(spacing) ')'], showOrient) ...
            repmat('local average orientions (red)', show_localOrient)] ... %repmat('NoLocalCoh', 1-show_localOrient) ...
            ...
           ['downsampling=' num2str(FactordsG) 'x' num2str(FactordsO) ... ', boxRadius=' num2str(boxRadius)  -- can be omitted with gradient
            ', coherency in box of radius ' num2str(boxRadius_LocalAlignment)], ...
           ['global coh=' num2str(ROIs_paths_struct(ll).overallCoherency) ...
            ', mean local coh=' num2str(ROIs_paths_struct(ll).meanLocalCoherency) ] } )
        
        %% zoom the figure to one ROI        
        %find coordinates of the ROI
        [yWhereROI, xWhereROI, ~] = find(maskFromROI);
        xlims = [min(xWhereROI), max(xWhereROI)] ;
        ylims = [min(yWhereROI), max(yWhereROI)] ;

        xlim(xlims)
        ylim(ylims)

        drawnow
        pause(0.1)


        % save files into one subfolder per ROI
        subfolderThisROI = [ output_fullPath(1:end-4) '_ROI' num2str(ll)];
        mkdir( subfolderThisROI )

        [~, filename_processingInfo, ~] = fileparts(filename_intermediate_output1);

        filename_print_singleROI = [subfolderThisROI  filesep  filename_processingInfo ...
            '_ROI' num2str(ll) '_' ... % only the subfolder and the added ROI number differ from the filename of full image plot
            repmat('Coherencyfield', showLocalCoherencyInsteadOfOrig) ...
            repmat('Orientations', showOrient) ...
            repmat('LocalCoh', show_localOrient) ... %repmat('NoLocalCoh', 1-show_localOrient) ...
            ... %NOTE: boxRadius_LocalAlignment is contained in filename_intermediate_output1
            repmat('CohOfROI',write_ROICoherency) ...
            repmat('OrientOfROI',show_ROIwiseOrient) ] ;


        print(Fig12, filename_print_singleROI,'-dpdf', '-r0', '-fillpage')
        print(Fig12, filename_print_singleROI,'-dpng', '-r0')
    end

    xlim('auto')
    ylim('auto')

    axis equal

    title({[repmat('Coherency field', showLocalCoherencyInsteadOfOrig) ...
            repmat('orig image', 1-showLocalCoherencyInsteadOfOrig) ...
        ' with ROIs outlined, ' ...
        repmat(['orientations (blue, shown every ' num2str(spacing) ')'], showOrient) ...
        repmat('local average orientions (red)', show_localOrient)] ... %repmat('NoLocalCoh', 1-show_localOrient) ...
        ...
       ['downsampling=' num2str(FactordsG) 'x' num2str(FactordsO) ... ', boxRadius=' num2str(boxRadius)  -- can be omitted with gradient
        ', coherency in box of radius ' num2str(boxRadius_LocalAlignment)]} )


    drawnow
    pause(0.1)


    filename_print = [filename_intermediate_output1(1:end-4) '_' ...
        repmat('Coherencyfield', showLocalCoherencyInsteadOfOrig) ...
        repmat('Orientations', showOrient) ...
        repmat('LocalCoh', show_localOrient) ... %repmat('NoLocalCoh', 1-show_localOrient) ...
        ... %NOTE: boxRadius_LocalAlignment is contained in filename_intermediate_output1
        repmat('CohOfROI',write_ROICoherency) ...
        repmat('OrientOfROI',show_ROIwiseOrient) ]


    print(Fig12, filename_print,'-dpdf', '-r0', '-fillpage')
    print(Fig12, filename_print,'-dpng', '-r0')


close all

end  %loop over different plot settings

end  %loop over images

end  %loop over folders within input folder, that is the loop over specimen ages
disp('ANALYSIS DONE')
datetime('now')



%% save collected coherencies
% do not overwrite files due to invested computation time 
filename_output_coherencies_collection = strrep(filename_coherencies_collection, '.txt', '.mat');
if ~exist(filename_output_coherencies_collection,'file')
    save( filename_output_coherencies_collection, ...
            'coherencies_collection_arr')
else
    versionNum = 1;
    while exist( strrep(filename_output_coherencies_collection, '.mat', ['_v' num2str(versionNum) '.mat']) , 'file')
        versionNum = versionNum+1;
    end
    save(        strrep(filename_output_coherencies_collection, '.mat', ['_v' num2str(versionNum) '.mat']) , ...
            'coherencies_collection_arr')
    warning(['Saving cell array  coherencies_collection_arr  to a new version #' num2str(versionNum)])
end


