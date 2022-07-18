
%%% The analytical Matlab code in this file was developed for a biomedical
%%% research project. The results of that project were published during 2022
%%% in the ACS Chemical Biology peer-reviewed scientific journal.

%%% All reported usage of this work should cite the below published
%%% article. See journal article for description of functions in this code.

%%% Article Title: "A real-time cellular thermal shift assay (RT-CETSA) to monitor target engagement"
%%% Article Authors: Sanchez, Tino; Ronzetti, Michael; Owens, Ashley; Antony, Maria; Voss, Ty; Wallgren, Eric; Talley, Daniel; Balakrishnan, Krishna; Leyes Porello, Sebastian; Rai, Ganesha; Marugan, Juan; Michael, Samuel; Baljinnyam, Bolormaa; Southall, Noel; Simeonov, Anton; Henderson, Mark J.
%%% Journal Title: ACS Chemical Biology
%%% Year Published: 2022


clear
% set directory that contains input tif files from time lapse imaging
InputDir = uigetdir

cd(InputDir)

% this sets size of region in image that will be measured
% for each individual sample within the grid/plate layout, modify if plate
% or instrument set-up is changed
SizeRegionPerWell_x = 40;
SizeRegionPerWell_y = 40;

% opens first image in the time series, allows image parameters to be determined
% file name can be changed if needed, ie: systematic naming convention
% changes for another instrument
Im_1 = imread('1.tif');

% determines if image is RGB format grayscale (3x 8-bit channels).
% extracts 1 channel for further analysis
if size(Im_1, 3) > 1
   Im_1 = Im_1(:,:,1); 
end

% find im values for display
Im_1_PixInt_Min = min(Im_1(:));
Im_1_PixInt_Max = max(Im_1(:));


% user selects points to define grid layout
h_f1 = figure, imshow(Im_1,[])

uiwait(msgbox('select center of upper left well with mouse then hit enter', 'modal'))
[UL_x,UL_y] = getpts(h_f1)

uiwait(msgbox('select center of upper right well with mouse then hit enter', 'modal'))
[UR_x, UR_y] = getpts(h_f1)

uiwait(msgbox('select center of lower right well with mouse then hit enter', 'modal'))
[LR_x, LR_y] = getpts(h_f1)

% plot points for confirmation
Im_1_points = false(size(Im_1));
Im_1_points(uint16(UL_y),uint16(UL_x))=1;
Im_1_points(uint16(UR_y),uint16(UR_x))=1;
Im_1_points(uint16(LR_y),uint16(LR_x))=1;

% figure, imshow(Im_1_points,[]), title('selected corner points')

PlateChoice_List = {'96 well','384 well','1536 well'};
[PlateChoice_Index,PlateChoice_tf] = listdlg('ListString',PlateChoice_List,'PromptString',{'Select your plate type.',''},'SelectionMode','single');

% NumTotalCol = 24;
% NumTotalRow = 16;

NumTotalCol_List = [12; 24; 48];

NumTotalCol_First = 1;
NumTotalCol_Last = NumTotalCol_List(PlateChoice_Index);

NumTotalRow_List = [8; 16; 32];

NumTotalRow_First = 1;
NumTotalRow_Last = NumTotalRow_List(PlateChoice_Index);


prompt = {'Enter first row coordinate selected:','Enter first column coordinate selected:','Enter last row coordinate selected:','Enter last column coordinate selected:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {num2str(NumTotalRow_First), num2str(NumTotalCol_First), num2str(NumTotalRow_Last), num2str(NumTotalCol_Last)};
answer_PlateCoord = inputdlg(prompt,dlgtitle,dims,definput);

NumTotalRow_First = str2num(cell2mat(answer_PlateCoord(1)));
NumTotalCol_First = str2num(cell2mat(answer_PlateCoord(2)));
NumTotalRow_Last = str2num(cell2mat(answer_PlateCoord(3)));
NumTotalCol_Last = str2num(cell2mat(answer_PlateCoord(4)));

Tif_File_List_1 = ls ('*.tif');

% number of timepoints based on number of image files in time series.
Num_TimePoints = size(Tif_File_List_1,1);
OutputNumbers_MeanInt_Uncorrected = uint16(zeros(((NumTotalCol_Last-NumTotalCol_First)*(NumTotalRow_Last-NumTotalRow_First)),Num_TimePoints));
OutputNumbers_MeanInt_SubtractBG = OutputNumbers_MeanInt_Uncorrected;

Mask_Outline_Single = false(SizeRegionPerWell_y+1,SizeRegionPerWell_x+1);
Mask_Outline_Single(:,1)=1;
Mask_Outline_Single(:,end)=1;
Mask_Outline_Single(1,:)=1;
Mask_Outline_Single(end,:)=1;

% loop for time points
for Count_Time = 1:Num_TimePoints
    % Count_Time = 1;
    
    Count_Time
    
    Output_PlateLayout_1 = uint16(zeros(NumTotalRow_List(PlateChoice_Index),NumTotalCol_List(PlateChoice_Index))); 
    
    
    myfilename = sprintf('%d.tif', Count_Time);
  
    Im_1 = imread(myfilename);
  
    % determines if image is RGB format grayscale (3x 8-bit channels).
    % extracts 1 channel for further analysis
    if size(Im_1,3) > 1
        Im_1 = Im_1(:,:,1);
    end

    % set output couter to zero for each grid position (biological sample)
    Count_Output = 0;
        
    % increase size of points for display of central well region position on image   
    Im_1_points_dil = imdilate(Im_1_points, strel('disk',1));
    % figure, imshow(Im_1_points_dil,[])

    % y scale is backward becuase it is an image, flips sign of slope
    Line_UL_UR_slope = (UL_y - UR_y) / (UL_x - UR_x);
    Line_UL_UR_degrees = atand(Line_UL_UR_slope);

    % corrects for any rotation of the image,
    % based on user selected sample positions
    Im_2 = imrotate(Im_1, Line_UL_UR_degrees);
    % figure, imshow(Im_2,[Im_1_PixInt_Min Im_1_PixInt_Max])

    % create temp image that can later be used 
    % for well center visualization, useful for confirmation or debugging
    test1 = Im_2;

    Im_2_points_dil = imrotate(Im_1_points_dil, Line_UL_UR_degrees);
    % figure, imshow(Im_2_points_dil,[])

    Im_2_points_dil_L = bwlabel(Im_2_points_dil);

    Im_2_points_props = regionprops(Im_2_points_dil_L, 'Centroid');



    % the below values should be left as floating point format... 
    % warning: loss of precision. should convert to integer at last step

    % Rotated_WellCenters_XY =
    % reshape([Im_2_points_props.Centroid],[size(Im_2_points_props,1),2]);
    % wrong format above
    Rotated_WellCenters_XY = reshape([Im_2_points_props.Centroid],[2,size(Im_2_points_props,1)])'; % correct matrix fomat for paired xy points


    % Rotated_WellCenters_PerPlate_Xmin = (Im_2_points_props(1).Centroid(1));
    % Rotated_WellCenters_PerPlate_Xmax = (Im_2_points_props(2).Centroid(1));

    Rotated_WellCenters_PerPlate_Xmin = min(Rotated_WellCenters_XY(:,1));
    Rotated_WellCenters_PerPlate_Xmax = max(Rotated_WellCenters_XY(:,1));

    Rotated_WellCenters_PerPlate_Xdist = Rotated_WellCenters_PerPlate_Xmax - Rotated_WellCenters_PerPlate_Xmin;

    Rotated_WellDistance_X = Rotated_WellCenters_PerPlate_Xdist / ((NumTotalCol_Last - NumTotalCol_First));


    % Rotated_WellCenters_PerPlate_Ymin = (Im_2_points_props(1).Centroid(2));
    % Rotated_WellCenters_PerPlate_Ymax = (Im_2_points_props(3).Centroid(2));

    Rotated_WellCenters_PerPlate_Ymin = min(Rotated_WellCenters_XY(:,2));
    Rotated_WellCenters_PerPlate_Ymax = max(Rotated_WellCenters_XY(:,2));


    Rotated_WellCenters_PerPlate_Ydist = Rotated_WellCenters_PerPlate_Ymax - Rotated_WellCenters_PerPlate_Ymin;

    Rotated_WellDistance_Y = Rotated_WellCenters_PerPlate_Ydist / ((NumTotalRow_Last - NumTotalRow_First));


    SizeRegionPerWell_x_HalfDist = (SizeRegionPerWell_x/2);
    SizeRegionPerWell_y_HalfDist = (SizeRegionPerWell_y/2);


    Mask_Outline_All = uint8(zeros(size(Im_2)));
    Count_Well_Per_Time = 0;

    CountRow = 0;
    CountCol = 0;

        for C1 = NumTotalRow_First:NumTotalRow_Last % sample grid row loop
            % C1 = 1
            % C1
            CountRow = CountRow + 1;

            for C2 = NumTotalCol_First:NumTotalCol_Last % sample grid column loop
              % C2 = 10
              % C2
              CountCol = CountCol + 1;

              Count_Well_Per_Time = Count_Well_Per_Time + 1;

              SingleWell_Center_X = uint16(((C2-1)*Rotated_WellDistance_X) + Rotated_WellCenters_PerPlate_Xmin);
              SingleWell_Center_Y = uint16(((C1-1)*Rotated_WellDistance_Y) + Rotated_WellCenters_PerPlate_Ymin);

              % test1 = Im_2;
              test1(uint16(SingleWell_Center_Y), uint16(SingleWell_Center_X)) = 40000;

              % 2018b
              % h_ROI = images.roi.Rectangle(gca,'Position',[(SingleWell_Center_X - SizeRegionPerWell_x_HalfDist),(SingleWell_Center_Y - SizeRegionPerWell_y_HalfDist),SizeRegionPerWell_x,SizeRegionPerWell_y],'StripeColor','r');
                
              % Extract single sample region and convert to 16 bit: makes data consistent across different
              % input types
              Im_2_SingleWell = Im_2(uint16((SingleWell_Center_Y - SizeRegionPerWell_y_HalfDist)):uint16((SingleWell_Center_Y + SizeRegionPerWell_y_HalfDist)), uint16((SingleWell_Center_X - SizeRegionPerWell_x_HalfDist)):uint16((SingleWell_Center_X + SizeRegionPerWell_x_HalfDist)));
              % figure, imshow(Im_2_SingleWell,[]), title(strcat('r', num2str(CountRow), ' ', 'c', num2str(CountCol)))

              % median filter to reduce outlier pixels prior to
              % thresholding
              Im_2_SingleWell_MedFilt = medfilt2(Im_2_SingleWell,[5 5], 'symmetric');
               % figure, imshow(Im_2_SingleWell_MedFilt,[])

              % threshold local region 
              [Thrsh,Thrsh_EM] = graythresh(Im_2_SingleWell_MedFilt);

              % select lowest 5% of pixels, then measure median to estimate
              % local background
              BG_x1 = sort(Im_2_SingleWell_MedFilt);
              BG_x2 = BG_x1(1:uint16(size(BG_x1,1)/20));
              BG_x3 = median(BG_x2(:));

              if Count_Well_Per_Time == 1 % set this on first time point only
                  BG_Well_perTime(Count_Well_Per_Time,1) = BG_x3;
              end
                
             % incremental counter for each sample position (row column combination) 
             Count_Output = Count_Output + 1;

             OutputNumbers(Count_Output,1) = C1; % well row coord
             OutputNumbers(Count_Output,2) = C2; % well column coord
            % OutputNumbers(Count_Output,3) = Mask_Area;
            % OutputNumbers(Count_Output,4) = Mask_MeanInt_Uncorrected;
            % OutputNumbers(Count_Output,5) = BG_MeanInt;
            % OutputNumbers(Count_Output,6) = Mask_MeanInt_SubtractBG;

            OutputNumbers_MeanInt_Uncorrected(Count_Output,Count_Time) = mean(Im_2_SingleWell(:));
            % OutputNumbers_MeanInt_SubtractBG(Count_Output,Count_Time) = Mask_MeanInt_SubtractBG;
            OutputNumbers_MeanInt_SubtractBG(Count_Output,Count_Time) = mean(Im_2_SingleWell(:))-BG_x3;

            OutputNumbers_BG(Count_Output,Count_Time)=BG_x3;

            % Mask_Outline_Single = bwperim(Mask_2);
            % Mask_Outline_Single = true(size(Im_2_SingleWell));

            Mask_Outline_Single_FullField = false(size(Im_2)); 
            Mask_Outline_Single_FullField(((SingleWell_Center_Y - SizeRegionPerWell_y_HalfDist)):uint16((SingleWell_Center_Y + SizeRegionPerWell_y_HalfDist)), uint16((SingleWell_Center_X - SizeRegionPerWell_x_HalfDist)):uint16((SingleWell_Center_X + SizeRegionPerWell_x_HalfDist))) = Mask_Outline_Single;

            Mask_Outline_All(Mask_Outline_Single_FullField == 1) = 200;
            % end % end conditional for threshold quality



            end % end C2 Column Loop
        end % C1 Row Loop
    
    % figure, imshow(Mask_Outline_All,[])
    
    % for visualizing/debugging sample center positions
    test1 = imdilate(test1, strel('disk',3));
    % imshow(test1,[min(Im_1(:)) max(Im_1(:))]), title(strcat('TimePoint ',num2str(Count_Time)))
    
    % for display
    ScaleInt =  double(max(Im_2(:))-median(BG_Well_perTime(:)))/256;
    Im_2_NoBG = Im_2 - double(median(BG_Well_perTime(:)));
    Im_2_NoBG_8bit = uint8(Im_2_NoBG/ScaleInt);
    % figure, imshow(Im_2_NoBG_8bit,[])
    
    % for display
    Im_2_RegionOverlay = uint8(Mask_Outline_All);
    Im_2_RegionOverlay(:,:,2)=Im_2_NoBG_8bit;
    Im_2_RegionOverlay(:,:,3) = 0;
      imshow(Im_2_RegionOverlay,[])
    
    
end  % end loop for time points



 % concatenate well coordinates with measurements
% % OutputList_SortedByTime = [OutputNumbers(:,1:2),OutputNumbers_MeanInt_SubtractBG];
% % OutputList_SortedByTime = [OutputNumbers(:,1:2),OutputNumbers_MeanInt_Uncorrected];
% % OutputList_SortedByTime = [OutputNumbers(:,1:2),OutputNumbers_BG];

% correct background based on first time point in the time series
OutputNumbers_SubtractSingleBG = double(OutputNumbers_MeanInt_Uncorrected) - double(median(BG_Well_perTime(:)));

% For output, combine sample row column coordinates with background subtracted results 
OutputList_SortedByTime = [OutputNumbers(:,1:2),OutputNumbers_SubtractSingleBG];

% save output as excel file in same directory as input files
Filename = 'MatlabResults_1.xlsx';
writematrix(OutputList_SortedByTime, Filename,'FileType','spreadsheet')


% 
%copy above text to analyze plate
%

%%% use below commented code to initially compare results for samples
%%% Plot curves for 2 wells
% Plot_Well_1_Coord_Row = 2
% Plot_Well_1_Coord_Col = 1
% 
% Plot_Well_2_Coord_Row = 5
% Plot_Well_2_Coord_Col = 3
% 
% Plot_Well_1_Idx = find(OutputList_SortedByTime(:,1) == Plot_Well_1_Coord_Row & OutputList_SortedByTime(:,2) == Plot_Well_1_Coord_Col);
% Plot_Well_2_Idx = find(OutputList_SortedByTime(:,1) == Plot_Well_2_Coord_Row & OutputList_SortedByTime(:,2) == Plot_Well_2_Coord_Col);
% 
% y1 = OutputList_SortedByTime(Plot_Well_1_Idx,3:end)';
% y2 = OutputList_SortedByTime(Plot_Well_2_Idx,3:end)';
% x1 = 1:size(y1,1)';
% 
% figure
% plot(x1,y1) 
% hold on 
% plot(x1,y2)
% legend(gca,'show');