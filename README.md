# segmentation-toolkit

The segmentation toolkit contains functions required for cell segmentation in drosophila embryo germband extension movies.

Cell segmentation is done using a [Watershed](https://www.mathworks.com/help/images/ref/watershed.html) transformation.

## Segmentation Flowchart (2D)
![](images/Segmentation_code_flowchart_2D.png)

## Segmentation steps (2D):
1. Create image file list.
    - Navigate in Matlab to the folder location where the movie images are stored. This will be used as the path in the ELSA1 function.
    - The ELSA1 function will bring up a finder window. You need to select the first z layer of the first time point, then the first z layer of the second time point.

    ```matlab
    path = cd
    [imageFileList] = ELSA_1loadImageList(path)
    ```

    - ELSA1 function is dependent on `loadImageNames_dualIndex`, `getFileStackNames_two`, `getFilenameBody`, and `getFilenameBody_two` functions.

2. Store list and other data in a structure.
    - Important data to save would be at a minimum the ImageFileList, Source, SecPerFrame. Other data can include start and end of GBE, pixel size in microns, processed movies, etc.
    - Save the data structure in the Source folder.

    ```matlab
    data.ImageFileList = imageFileList
    data.Source = path
    data.SecPerFrame = 1
    save('data','data')
    ```

    - If the movie has been segmented and analyzed previously but you need to change the ImageFileList and Source, you can use the function `rewriteImageFileListAndSource`.
    - The 'oldString' and 'newString' inputs refer to portions of the path that need to be changed, the whole path does not need to be inputted.
    ```matlab
    [dataNew] = rewriteImageFileListAndSource(data,'oldString','newString')
    ```

3. Initialize seeds for the first movie frame.
    - The ELSA2 function will by default pick the first frame to initialize seeds, this can be changed in the code inputs.
    - The function will bring up an interactive image of seeds and mask which can be edited using the mouse and keyboard.
    ```matlab
    ELSA_2initializeSeedsAndMask(data)
    ```
    - ELSA2 function is dependent on `im2seedsAndMask` (local), `modifySeedsAndMask_ForInitialization` (local), and `filterImage3DpaddedEdges` functions.
4. Run segmentation and check the results every so often. Add seeds when cells are sufficiently in the frame.
    - The ELSA3 function gives options for the time period to segment (given in frames, i.e. [1:20]), as well as the check interval (number of frames before checking, i.e. 5). These will change depending on the movie to be segmented.
    - The function will bring up an interactive image of seeds and mask which can be edited. If there is a mistake, there is an "Undo" option which will revert the seeds of that frame to their original configuration.

    ```matlab
    ELSA_3segmentCellsV5(data,5,[1:20])
    ```

    - ELSA3 function is dependent on `wsSegmentSingleImageV5` (local), `modfiySeedsAndMaskV5` (local), `im2seeds` (local), and `filterImage3DpaddedEdges` functions. 
5. Run cell tracking.
    - The ELSA4 function will take the labeled cell images and track them over time so cell IDs are unique for the whole movie.
    - This function assumes using only one z-layer (2D segmentation).

    ```matlab
    ELSA_42trackCellsInT_SingleZ(data)
    ```

    - There is a function called that is written in C++ (`mexLap`). In order for it to run properly, the correct library must be used. For example, if running on a mac the extension `mexLap.mexmaci64` must be available on the computer. If the library is not available, you can compile it in Matlab using the mex function. For example, in Matlab the commands might look like: 

    ```matlab
    mex -setup C++
    mex mexLap.cpp
    ```


    - ELSA4 function is dependent on `trackArea2`, `linkFrame2Frame_ver2` (local), `conflictResolution` (local), `lap`, and `mexLap` functions.


6. Run node extraction.
    - The ELSA5 function gives the option for the time period to segment (given in frames, i.e. [1:20]). If this is not inputted, it will use the min and max available. It is important to note that if you segmented less than the full amount of available frames, you may run into an error if the time vector is higher.

    ```matlab
    ELSA_52extractNodes_V2(data,[1:20])
    ```

    - ELSA5 function is dependent on `nodeAnalysis`, `detectNodes_ver2`, `DistanceMatrix` (local), and `linkNodesAndcellsV2` functions.

7. Collect interface based data using the gridAnalysis function.
    - The gridAnalysis function gives the option for the first and last frame (i.e. [1 20]). If this is not inputted, it will use the min and max available. It is important to note that if you segmented less than the full amount of available frames, you may run into an error if the time vector is higher.
    - The gridAnalysis function gives the option for an angle vector. If you want all angles, this can be inputted as empty (i.e. []).

    ```matlab
    [matRes,trackingMatrixZT,trackingVectorLI] = gridAnalysis_interfaceOrientationNewZT6(data.Source,[], data.ImageFileList,[1 20])
    ```

    - gridAnalysis function is dependent on `thresholdVector` function.


8. Collect cell centroid and cell geometry (area and perimeter) data.
    - The centroidArray and GeometryData functions give the option for a movie number in case the data structure stores information for multiple movies. If only one movie is stored in the data structure, the input should just be 1.
    ```matlab
    centroidArray(data,1)
    ```
    ```matlab
    GeometryData(data,1)
    ```


9. Collect other data for additional analysis (i.e. intensity data).
    - This code is not included in the segmentation-toolkit.


10. Check results from segmentation and tracking as needed.
    - Results from segmentation can be checked retrospectively if needed (i.e. if downstream analysis is giving strange results).
    - The viewresults function will allow for checking the raw image, the tracking labels, the segmentation image, the segmentation lines overlaid on raw image, or the seeds and mask.
    - This function gives the option for frame delay (time interval between frames played back) and the time vector (specific frames to analyze, i.e. [1:20]).
    ```matlab
    viewresults(data,'raw',1,[1:20])
    ```
    - viewresults function is dependent on `imoverlay` (local) function.

## Segmentation Flowchart (3D)
![](images/Segmentation_code_flowchart_3D.png)

## Segmentation Steps (3D)

1. Create image file list.
    - Navigate in Matlab to the folder location where the movie images are stored. This will be used as the path in the ELSA1 function.
    - The ELSA1 function will bring up a finder window. You need to select the first z layer of the first time point, then the first z layer of the second time point.

    ```matlab
    path = cd
    [imageFileList] = ELSA_1loadImageList(path)
    ```

    - ELSA1 function is dependent on `loadImageNames_dualIndex`, `getFileStackNames_two`, `getFilenameBody`, and `getFilenameBody_two` functions.

2. Store list and other data in a structure.
    - Important data to save would be at a minimum the ImageFileList, Source, SecPerFrame, Nlayers. Other data can include start and end of GBE, pixel size in microns, processed movies, etc.
    - Save the data structure in the Source folder.

    ```matlab
    data.ImageFileList = imageFileList
    data.Source = path
    data.SecPerFrame = 1
    data.Nlayers = 20
    save('data','data')
    ```

    - If the movie has been segmented and analyzed previously but you need to change the ImageFileList and Source, you can use the function `rewriteImageFileListAndSource`.
    - The 'oldString' and 'newString' inputs refer to portions of the path that need to be changed, the whole path does not need to be inputted.
    ```matlab
    [dataNew] = rewriteImageFileListAndSource(data,'oldString','newString')
    ```

3. Create empty matrices for segmentation.
    - Downstream 3D segmentation functions pull empty matrices for seeds, mask, ImageSegment, and ImageBWlabel, which are then filled in subsequent steps. Therefore, these arrays must be initialized at the beginning.
    - This function gives the option for the time period to initialize arrays (given in frames, i.e. [1:20]), as well as the movie number (if there is only one movie in the data structure, this will just be 1).

    ```matlab
    ELSA_initializeEmptySegData_LS(data,1,[1:20])
    ```
4. Initialize seeds for the first movie frame.
    - This function has an input to specify the time frame to initialize seeds (if the first frame looks good, then this could be 1). It also has an input to specify the z layer to initialize the seeds (if the top z layer looks good, then this could be 1).
    - This function has an input to specify the seedfilter to use. This is a gaussian filter used for filtering the image. 
    - The function will bring up an interactive image of seeds and mask which can be edited using the mouse and keyboard.

    ```matlab
    ELSA_initialize3Dseeds_SingleLayer_LS(data,1,1,seedfilter)
    ```
    - This function is dependent on `modifySeedsAndMask_LS` (local), `im2seeds` (local) and `filterImage3DpaddedEdges` functions.
5. Run cell segmentation and tracking in the time dimension.
    - This function will segment and track images over time so cell IDs are unique for the whole movie.
    - This function uses only one z-layer to track through time. The z-layer inputted should be the same layer that was used in the initialization of seeds.
    - The function give the option for a movie number in case the data structure stores information for multiple movies. If only one movie is stored in the data structure, the input should just be 1.
    - The function gives options for the time period to segment (given in frames, i.e. [1:20]), as well as the check interval (number of frames before checking, i.e. 5). These will change depending on the movie to be segmented.

    ```matlab
    ELSA_segmentCellsInTime_LS(data,1,5,[1:20],1)
    ```
    - This function is dependent on `labelTracker` (local), `wsSegmentSingleImageV5` (local), `modifySeedsAndMask_LS` (local), `im2seeds` (local), and `filterImage3DpaddedEdges` functions.

6. Propagate segmentation in the z-dimension.
    - This function will take the segmented frames that have been tracked over time and propagate the segmentation through all z-layers.
    - This function asks for a starting z-layer to be inputted. This should be the same layer that was used in the initialization of seeds and used for segmenting in time.
    - The function gives the option for a movie number in case the data structure stores information for multiple movies. If only one movie is stored in the data structure, the input should just be 1.
    - The function gives options for the time period to segment (given in frames, i.e. [1:20]).

    ```matlab
    ELSA_segmentCells_PropagateInZ_LSV2(data,1,2,[1:20])
    ```

    - This function is dependent on `labelTracker` (local) and `wsSegmentSingleImageV5` (local) functions.
7. Run node extraction. 
    - This function gives the option for the time period to segment (given in frames, i.e. [1:20]). If this is not inputted, it will use the min and max available. It is important to note that if you segmented less than the full amount of available frames, you may run into an error if the time vector is higher.

    ```matlab
    ELSA_52extractNodes_V3(data,[1:20])
    ```

    - ELSA5 function is dependent on `nodeAnalysis`, `detectNodes_ver2`, `DistanceMatrix` (local), and `linkNodesAndcellsV2` functions.

8. Collect interface based data using the gridAnalysis function.
    - The gridAnalysis function gives the option for the first and last frame (i.e. [1 20]). If this is not inputted, it will use the min and max available. It is important to note that if you segmented less than the full amount of available frames, you may run into an error if the time vector is higher.
    - The gridAnalysis function gives the option for an angle vector. If you want all angles, this can be inputted as empty (i.e. []).

    ```matlab
    [matRes,trackingMatrixZT,trackingVectorLI] = gridAnalysis_interfaceOrientationNewZT6(data.Source,[], data.ImageFileList,[1 20])
    ```

    - gridAnalysis function is dependent on `thresholdVector` function.

9. Collect cell centroid and cell geometry (area and perimeter) data.
    - The centroidArray3D and GeometryData functions give the option for a movie number in case the data structure stores information for multiple movies. If only one movie is stored in the data structure, the input should just be 1.
    - The GeometryData function requires the number of frames to be analyzed as an input (i.e. 20 to analyze from frame 1-20).
    ```matlab
    centroidArray3D(data,1)
    ```
    ```matlab
    GeometryData_LS(data,1,20)
    ```

10. Collect other data for additional analysis (i.e. intensity data).
    - This code is not included in the segmentation-toolkit.


11. Check results from segmentation and tracking as needed.
    - Results from segmentation can be checked retrospectively if needed (i.e. if downstream analysis is giving strange results).
    - The viewresults3D function will allow for checking the raw image, the tracking labels, the segmentation image, the segmentation lines overlaid on raw image, or the seeds and mask.
    - This function gives the option for frame delay (time interval between frames played back) and the time vector (specific frames to analyze, i.e. [1:20]).
    - This function gives the option for z-layers to analyze (i.e. [1:5] if you want to look at z layers 1-5).
    ```matlab
    viewresults3D(data,'raw',1,[1:20],[1])
    ```
    - viewresults3D function is dependent on `imoverlay` (local) function.