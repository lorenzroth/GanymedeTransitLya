close all;
clear all;

addpath("../src")
addpath("../utility_files") 
addpath("../ganymede_observations");
run("loadParametersAnalysis.m") % load parameters


panel_figure = figure("Units","inches","Position",[1,1,10,8]);
counter = 1;

for observation_ID= [3010,4010]

    which_observation = observation_ID;

    if which_observation == 3010
    observation = observation3010;

    elseif which_observation == 4010
        observation = observation4010;
    
    elseif which_observation == 3011
        observation = observation3011;
    
    elseif which_observation == 4011
        observation = observation4011;
    
    else
        error("uknown obseravtion")
    end

    observations_dir   = "../ganymede_observations";

    % extract fits data
    GanymedeImage = FitsImageObject;
    GanymedeImage.read_image(fullfile(observations_dir,observation.file));

    %%  PLOT IMAGE
    [rows,cols]              = size(GanymedeImage.image);     % find image size
    x_pixel_range_full_image = 1:cols;                        % full image axis range  x
    y_pixel_range_full_image = 1:rows;                        % full images axis range y
    max_intensity            = max(max(GanymedeImage.image)); % define max intensity
    %cscale                   = [0,max_intensity];
    cscale                   = [0,20];

    % dark counts should remain for the analysis, but you don't want
    % to see them in the picture

    ax   = subplot(2,2,counter);
    imagesc(ax,GanymedeImage.image,cscale);
    hold on
    axis square;
    set(gca,'ydir','normal') ;  % y-axis should go from bottom to top (vertical)
    ax.XLim = [1,cols];
    ax.YLim = [1,rows];
    ax.YLabel.String = "y pixel";
    ax.YTick = 100:200:1000;
    ax.XTick = 100:200:1000;
    ax.Title.String  = observation.name;
    ax.XMinorTick = "on";
    ax.TickLength= [0.010,0.010]
    ax.LineWidth = 2;
    ax.Box = "on";
    ax.BoxStyle ="full";
    color_bar_instance = formal_colorbar(colorbar(ax));
    color_bar_instance.Label.String = "Detector Counts";
    color_bar_instance.Ticks = 0:5:40;
    %clim([0 40])

    ax_under = doubleXAxis(ax);
    hold on
    axis square;
    [x_wavelength,~] = GanymedeImage.obtain_axis_conversions();
    ax_under.XLim = [x_wavelength(1),x_wavelength(end)];
    ax_under.XTick = 1100:100:1700;
    ax_under.XLabel.String = "x pixel / Wavelength (A)";

    %% PREPARE PLOT FOR SUBIMAGE

    % define subimage outer edge
    x_pmin = observation.x_pcenter-95 ; % pixel
    x_pmax = observation.x_pcenter+95; % pixel
    y_pmin = observation.y_pcenter-95; % pixel
    y_pmax = observation.y_pcenter+95; % pixel

    % plot borders
    x_square = [x_pmin,x_pmin,x_pmax,x_pmax];
    y_square = [y_pmin,y_pmax,y_pmin,y_pmax];
%     scatter(ax,x_square,y_square,'green','filled',LineWidth=10)

    %% RESIZE IMAGE
    [ganymede_centred_subimage,sigma_matrix_ganymede_centred_subimage] = GanymedeImage.resize_image(x_pmin,x_pmax,y_pmin,y_pmax);
    x_pixel_range_ganymede_centred_subimage     = x_pixel_range_full_image(x_pixel_range_full_image<=x_pmax & x_pixel_range_full_image>= x_pmin);
    y_pixel_range_ganymede_centred_subimage     = y_pixel_range_full_image(y_pixel_range_full_image<=y_pmax & y_pixel_range_full_image>= y_pmin);


    %% CENTER DEFINITION

    x_pixel_center       = observation.x_pcenter ;   
    y_pixel_center       = observation.y_pcenter ;   
    box_radial_extension = 1.5;      % box around the ganymede to be eliminated
    % from the fit expressed is ganymedes radii

    x_index_center_ganymede_centred_subimage = find(x_pixel_range_ganymede_centred_subimage == x_pixel_center); % index of the center in the x_range_sub
    y_index_center_ganymede_centred_subimage = find(y_pixel_range_ganymede_centred_subimage == y_pixel_center); % index of the center in the y_range_sub

    %% CONVERSION FROM COUNTS TO REYLIGHTS

    exposition_time =  GanymedeImage.find_key("TEXPTIME"); %s

    filter_data   = dlmread("HST_STIS_FUV.25MAMA_G140L.dat");
    wavelength    = filter_data(:,1);
    throughput    = filter_data(:,2);
    throughput_Ly = interp1(wavelength,throughput,1216);
    A_eff         = globalParameters.A_mirror*throughput_Ly;
    Omega         = globalParameters.mx*globalParameters.my*(2*pi/3600/360)^2;

    count2KRayleight = 4*pi/10^6/(exposition_time*Omega*A_eff)*10^-3;
    ganymede_centred_subimage_reyleights = ganymede_centred_subimage*count2KRayleight;
    sigma_matrix_ganymede_centred_subimage_reyleights = sigma_matrix_ganymede_centred_subimage*count2KRayleight;

    
    xcorner = [x_pixel_range_ganymede_centred_subimage(1),x_pixel_range_ganymede_centred_subimage(end)];
    ycorner  = [y_pixel_range_ganymede_centred_subimage(1),y_pixel_range_ganymede_centred_subimage(end)];

    max_intensity            = max(max(ganymede_centred_subimage_reyleights)); % define max intensity
    %cscale                   = [0,max_intensity];
    %cscale                   = [0,0.5]; for 1356 aurora test
    cscale                   = [4,10];

    ax   = subplot(2,2,counter+2);
    imagesc(ax,xcorner,ycorner,ganymede_centred_subimage_reyleights,cscale);
    hold on
    axis square;
    set(gca,'ydir','normal') ;  % y-axis should go from bottom to top (vertical)
    ax.XLim = [x_pixel_range_ganymede_centred_subimage(1),x_pixel_range_ganymede_centred_subimage(end)];
    ax.YLim = [y_pixel_range_ganymede_centred_subimage(1),y_pixel_range_ganymede_centred_subimage(end)];
    ax.XLabel.String = "x pixel"
    ax.YLabel.String = "y pixel";
    ax.YTick = 100:25:1000;
    ax.XTick = 100:25:1000;
    ax.Title.String  = observation.name;
    ax.XMinorTick = "on";
    ax.YMinorTick = "on";
    ax.TickLength= [0.010,0.010]
    ax.LineWidth = 2;
    ax.Box = "on";
    ax.BoxStyle ="full";
    color_bar_instance = formal_colorbar(colorbar());
    color_bar_instance.Label.String = "Brightness [kR]";
    color_bar_instance.Ticks = 0:2:30;
    %clim([0 25]);

    draw_circle([observation.x_pcenter,observation.y_pcenter],globalParameters.diameter_ganymede/2,ax)
    
    counter = counter+1;

    %fontsize(gca,scale=1.5)    

    
end