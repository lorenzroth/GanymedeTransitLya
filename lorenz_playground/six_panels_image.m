close all;
clear all;
addpath("../src")
addpath("../utility_files")

% THIS SCRIPT IS ALL AUTOMATIC AND IT RUNS THE 9 PANELS IMAGE FOR BOTH THE 
% FIGURES. YOU CAN SET THE PARAMETERS FOR ALL THE FIGURES AND THEN IT IS

panel_figure = figure("Position",[1,1,600,800]);
counter             = 1;
run("loadParametersAnalysis.m") % load parameters

for ID=[3011,4011]
   which_observation = ID;

if which_observation == 3010
    observation = observation3010;

elseif which_observation == 4010
    observation = observation4010;

elseif which_observation == 3011
    observation = observation3011;

elseif which_observation == 4011
    observation = observation4011;

else
    error("uknown observation")
end


if which_observation == 3010
    observation = observation3010;

elseif which_observation == 4010
    observation = observation4010;

elseif which_observation == 3011
    observation = observation3011;

elseif which_observation == 4011
    observation = observation4011;

else
    error("uknown observation")
end


observations_dir   = "../ganymede_observations";
filename           = observation.file;
output_dir         = "../images";
output_dir_models  = "../images/models";

thickness_annulus = 3; % in pixels 

% extract fits data
GanymedeImage = FitsImageObject;
GanymedeImage.read_image(fullfile(observations_dir,filename));

%%  PLOT IMAGE

[rows,cols]              = size(GanymedeImage.image);     % find image size
x_pixel_range_full_image = 1:cols;                        % full image axis range  x
y_pixel_range_full_image = 1:rows;                        % full images axis range y
max_intensity            = max(max(GanymedeImage.image)); % define max intensity
cscale                   = [0,max_intensity];

% dark counts should remain for the analysis, but you don't want
% to see them in the picture


% define subimage outer edge
x_pmin = observation.x_pcenter-95 ; % pixel
x_pmax = observation.x_pcenter+95; % pixel
y_pmin = observation.y_pcenter-95; % pixel
y_pmax = observation.y_pcenter+95; % pixel

% plot borders
x_square = [x_pmin,x_pmin,x_pmax,x_pmax];
y_square = [y_pmin,y_pmax,y_pmin,y_pmax];


%% RESIZE IMAGE
[ganymede_centred_subimage,sigma_matrix_ganymede_centred_subimage] = GanymedeImage.resize_image(x_pmin,x_pmax,y_pmin,y_pmax);
x_pixel_range_ganymede_centred_subimage     = x_pixel_range_full_image(x_pixel_range_full_image<=x_pmax & x_pixel_range_full_image>= x_pmin);
y_pixel_range_ganymede_centred_subimage     = y_pixel_range_full_image(y_pixel_range_full_image<=y_pmax & y_pixel_range_full_image>= y_pmin);


%% CENTER DEFINITION

% floor(diameter/2) and center by eye (referenced to the full image)
x_pixel_center       = observation.x_pcenter ;      % define center of the image by eye on the big image
y_pixel_center       = observation.y_pcenter ;    % define center of the image by eye on the big image
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


%% FIT THE MODEL
% assumed ganymede brightness
globalParameters.ganymede_assumed_brightness = 1.3; % KReyleights

% decide max order for background polynomial surface
poly_order = observation.poly_order;

syms n0
syms a [poly_order+1,poly_order+1] 
syms x y

polynomial_surface = 0;
for ii=0:poly_order
    for jj =0:poly_order-ii
        polynomial_surface = polynomial_surface + a(ii+1,jj+1)*x^ii*y^jj;
    end
end

zc        = observation.zc;
R_ganymede_pixel  = globalParameters.diameter_ganymede/2;

 % H corona absorbtion model outside ganymede disk
Nh  = n0*globalParameters.RGanymedeCm *(R_ganymede_pixel./sqrt(x.^2+y.^2)).*pi;
tau = zc*Nh;
T   = exp(-tau);

% H corona emission model outside ganymede disk
g_factor  = observation.g_factor; %[kR/cm2]
I_H_corona_emission = Nh*g_factor;

% Define point spread function
PSF = imread("fuvmama_1216_00.fits"); % point spread function
[rows_psf,cols_psf] = size(PSF);

% Expand axis definition to eliminate noise after convolution
right_expansion =  max(x_pixel_range_ganymede_centred_subimage)+1:max(x_pixel_range_ganymede_centred_subimage)+cols_psf*2;   % 2 is just a sefty margin
left_expansion  =  (min(x_pixel_range_ganymede_centred_subimage)-cols_psf*2):min(x_pixel_range_ganymede_centred_subimage)-1; % 2 is just a sefty margin
up_expansion    =  max(y_pixel_range_ganymede_centred_subimage)+1:max(y_pixel_range_ganymede_centred_subimage)+rows_psf*2;   % 2 is just a sefty margin
down_expansion  =  (min(y_pixel_range_ganymede_centred_subimage)-cols_psf*2):min(y_pixel_range_ganymede_centred_subimage)-1; % 2 is just a sefty margin

x_expanded_axis = [left_expansion,x_pixel_range_ganymede_centred_subimage,right_expansion];
y_expanded_axis = [down_expansion,y_pixel_range_ganymede_centred_subimage,up_expansion];
x_original_mask = logical([zeros(size(left_expansion)),ones(size(x_pixel_range_ganymede_centred_subimage)),zeros(size(right_expansion))]);
y_original_mask = logical([zeros(size(down_expansion)),ones(size(y_pixel_range_ganymede_centred_subimage)),zeros(size(up_expansion))]);


% create grid for bightness extraction
[rows,cols] = size(ganymede_centred_subimage_reyleights);

%create pixel mask for ganymede_centred image with (0,0) at ganymede center
[X_gird_ganymede_centred_subimage,Y_gird_ganymede_centred_subimage] = meshgrid(x_pixel_range_ganymede_centred_subimage- x_pixel_center,y_pixel_range_ganymede_centred_subimage-y_pixel_center);

%create pixel mask for expanded image
[X_grid_model,Y_grid_model] = meshgrid(x_expanded_axis- x_pixel_center ,y_expanded_axis- y_pixel_center );

brightness_mask                          = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2) < R_ganymede_pixel/2; % only for finding a brightness value take the brightness in half the radius of the ganymede
ganymede_mask_for_ganymede_centred_subimage = sqrt(X_gird_ganymede_centred_subimage.^2 + Y_gird_ganymede_centred_subimage.^2) < R_ganymede_pixel;   % mask covering the ganymede disk in ganymede_centred image
ganymede_mask_for_model_image            = sqrt(X_grid_model.^2 + Y_grid_model.^2) < R_ganymede_pixel;   % this the mask covering the ganymede in the expanded image

mean_brightness_ganymede = mean(mean(ganymede_centred_subimage_reyleights(brightness_mask)));
IPMandGEO                = mean_brightness_ganymede - globalParameters.ganymede_assumed_brightness;
IPMandGEO_counts         = IPMandGEO/count2KRayleight;

% adjust sigma values
sigma_matrix_ganymede_centred_subimage_reyleight   = sqrt(sigma_matrix_ganymede_centred_subimage.^2+IPMandGEO_counts)*count2KRayleight;

mean_brightness_ganymede = mean_brightness_ganymede - IPMandGEO;
ganymede_centred_subimage_reyleights  = ganymede_centred_subimage_reyleights - IPMandGEO;

% obtain polynomial fit of the image
% the x and y axis must be centered with zero at the ganymede center.
% Otherwise the transmission model won't work
[x_fit,y_fit,z_fit,weights_fit] = prepareSurfaceData(X_gird_ganymede_centred_subimage,Y_gird_ganymede_centred_subimage,ganymede_centred_subimage_reyleights,1./sigma_matrix_ganymede_centred_subimage_reyleight.^2);

exclude_set = ~excludedata(x_fit,y_fit,'box',[-box_radial_extension*R_ganymede_pixel ,...
                                              +box_radial_extension*R_ganymede_pixel ,...
                                              -box_radial_extension*R_ganymede_pixel ,...
                                              + box_radial_extension*R_ganymede_pixel]);
% backgroumd surface model
total_model_sym     = polynomial_surface.*T + I_H_corona_emission  ; % background surface plus the Transmission plus corona emission
total_model         = matlabFunction(total_model_sym);
fit_parameters      = string(symvar(total_model_sym));

% eliminate variables that are not to be fit in the model
eliminate_xyn0_mask = string(symvar(total_model_sym)) ~= 'x' & ...
                      string(symvar(total_model_sym)) ~= 'y' &  ...
                      string(symvar(total_model_sym)) ~= 'n0';

fit_parameters      = cellstr(fit_parameters(eliminate_xyn0_mask));

myfittype = fittype(total_model,...
                     'dependent',{'z'},'independent',{'x','y'},...
                     'coefficients',fit_parameters ,'problem','n0');
options = fitoptions(myfittype);

coefficients_upperbound = +10000;
coefficients_lowerbound = -10000;

options.Upper      = ones(1,length(fit_parameters))*coefficients_upperbound;
options.Lower      = ones(1,length(fit_parameters))*coefficients_lowerbound;
mid_point_bound    = (coefficients_upperbound+coefficients_lowerbound)/2;
options.Exclude    = exclude_set;
options.StartPoint = zeros(1,length(fit_parameters))*mid_point_bound;
% options.Weights    = weights_fit;


shaded_error_alphas = 0.4;
observation_color   = "r";
error_bar_color     = [.7 .7 .7]; % gray

line_style_best_fit   = "-";
line_style_zero_nzero = "--";
line_style_model_only = "-.";
line_style_array = [line_style_zero_nzero,line_style_best_fit];
counter2 = 1;
for n0_value=[0,observation.best_n0]
   
    model_fit = fit([x_fit,y_fit],z_fit,myfittype,options,"problem",n0_value);
    % create model image from fit
    bkg_only_model_image = double(model_fit(X_grid_model,Y_grid_model));
   
  
    % assume constant brightness at disk as constant
    row_model_with_ganymede = bkg_only_model_image;
    row_model_with_ganymede(ganymede_mask_for_model_image)  = mean_brightness_ganymede;  
    
    final_model_image_PSF   = conv2(row_model_with_ganymede,PSF,'same');% apply PSF
    bkg_only_model_image    = bkg_only_model_image(y_original_mask,x_original_mask);
    final_model_image_PSF   = final_model_image_PSF(y_original_mask,x_original_mask);
    rsw_model_image         = row_model_with_ganymede(y_original_mask,x_original_mask); % image size before the PSF convolution is applied
    
    
    %% GENERATE RADIAL PLOTS
    
    if counter2 ==1 

      [ganymede_centred_subimage_vertical_box_plot,~]   = generate_vertical_rectangle_plot(ganymede_centred_subimage_reyleights,...
                                                                                            x_index_center_ganymede_centred_subimage,...
                                                                                            y_index_center_ganymede_centred_subimage,...
                                                                                            globalParameters.diameter_ganymede);
  
      [ganymede_centred_subimage_horizontal_box_plot,~] = generate_horizontal_rectangle_plot(ganymede_centred_subimage_reyleights,...
                                                                                             x_index_center_ganymede_centred_subimage,...
                                                                                             y_index_center_ganymede_centred_subimage,...
                                                                                             globalParameters.diameter_ganymede);
  
      [ganymede_centred_subimage_radial_plot,~]         = generate_radial_plot(ganymede_centred_subimage_reyleights,...
                                                                               x_index_center_ganymede_centred_subimage,...
                                                                               y_index_center_ganymede_centred_subimage,...
                                                                               R_ganymede_pixel*2.5,...
                                                                               globalParameters.thickness_annulus);

      [bkg_vertical_box_plot,~]   = generate_vertical_rectangle_plot(bkg_only_model_image, ...
                                                                     x_index_center_ganymede_centred_subimage,...
                                                                     y_index_center_ganymede_centred_subimage,...
                                                                     globalParameters.diameter_ganymede);
  
      [bkg_horizontal_box_plot,~] = generate_horizontal_rectangle_plot(bkg_only_model_image,...
                                                                       x_index_center_ganymede_centred_subimage,...
                                                                       y_index_center_ganymede_centred_subimage,...
                                                                       globalParameters.diameter_ganymede);
  
        
      [bkg_radial_plot,~]         = generate_radial_plot(bkg_only_model_image,...
                                                         x_index_center_ganymede_centred_subimage,...
                                                         y_index_center_ganymede_centred_subimage,...
                                                         R_ganymede_pixel*2.5,globalParameters.thickness_annulus );

    end
    
    [model_image_vertical_box_plot,vertical_range,error_bar_vertical_box_plot]   = generate_vertical_rectangle_plot(final_model_image_PSF,...
                                                                                                                     x_index_center_ganymede_centred_subimage,...
                                                                                                                     y_index_center_ganymede_centred_subimage,...
                                                                                                                     globalParameters.diameter_ganymede,...
                                                                                                                     sigma_matrix_ganymede_centred_subimage_reyleight );

    [model_image_horizontal_box_plot,horizontal_range,error_bar_horizontal_box_plot] = generate_horizontal_rectangle_plot(final_model_image_PSF,...
                                                                                                                           x_index_center_ganymede_centred_subimage,...
                                                                                                                           y_index_center_ganymede_centred_subimage,...
                                                                                                                           globalParameters.diameter_ganymede,...
                                                                                                                           sigma_matrix_ganymede_centred_subimage_reyleight );
    [model_image_radial_plot,radial_range,error_bar_radial_plot]  = generate_radial_plot(final_model_image_PSF,...
                                                                                         x_index_center_ganymede_centred_subimage,...
                                                                                         y_index_center_ganymede_centred_subimage,...
                                                                                         R_ganymede_pixel*2.5,...
                                                                                         globalParameters.thickness_annulus,...
                                                                                         sigma_matrix_ganymede_centred_subimage_reyleight );
    


     % exclude disk points mask 
    vert_mask_nodik = vertical_range<= -R_ganymede_pixel | vertical_range >= R_ganymede_pixel ;
    hor_mask_nodik =  horizontal_range<= -R_ganymede_pixel | horizontal_range >= R_ganymede_pixel ;
    rad_mask_nodisk = radial_range >= R_ganymede_pixel;

    ganymede_centred_subimage_vertical_box_plot_nodisk   = ganymede_centred_subimage_vertical_box_plot(vert_mask_nodik);
    ganymede_centred_subimage_horizontal_box_plot_nodisk = ganymede_centred_subimage_horizontal_box_plot(hor_mask_nodik);
    ganymede_centred_subimage_radial_plot_nodisk         = ganymede_centred_subimage_radial_plot(rad_mask_nodisk);
    
    model_image_vertical_box_plot_nodisk   = model_image_vertical_box_plot(vert_mask_nodik);
    model_image_horizontal_box_plot_nodisk = model_image_horizontal_box_plot(hor_mask_nodik);
    model_image_radial_plot_nodisk         = model_image_radial_plot(rad_mask_nodisk);
    
    error_bar_vertical_box_plot_nodisk   = error_bar_vertical_box_plot(vert_mask_nodik);
    error_bar_horizontal_box_plot_nodisk = error_bar_horizontal_box_plot(hor_mask_nodik);
    error_bar_radial_plot_nodisk     = error_bar_radial_plot(rad_mask_nodisk);

    vertical_range_nodisk   = vertical_range(vert_mask_nodik);
    horizontal_range_nodisk = horizontal_range(hor_mask_nodik);
    radial_range_nodisk     = radial_range(rad_mask_nodisk);
    
    % Chi2 per image : radial, vertical and horizontal
    chi2_rad     = sum((model_image_radial_plot_nodisk         - ganymede_centred_subimage_radial_plot_nodisk        ).^2./error_bar_radial_plot_nodisk.^2  )/numel(model_image_radial_plot_nodisk);
    chi2_vert    = sum((model_image_vertical_box_plot_nodisk   - ganymede_centred_subimage_vertical_box_plot_nodisk  ).^2./error_bar_vertical_box_plot_nodisk.^2  )/numel(model_image_vertical_box_plot_nodisk  );
    chi2_hor     = sum((model_image_horizontal_box_plot_nodisk - ganymede_centred_subimage_horizontal_box_plot_nodisk).^2./error_bar_horizontal_box_plot_nodisk.^2 )/numel(model_image_horizontal_box_plot_nodisk );
    
    
    
    if counter2 == 1
       boxing_ax1 = formal_axes(subplot(3,3,counter));
       hold on
       plot(boxing_ax1,vertical_range/R_ganymede_pixel,ganymede_centred_subimage_vertical_box_plot,"DisplayName","STIS observation")
       plot(boxing_ax1,vertical_range/R_ganymede_pixel,bkg_vertical_box_plot,"color","k","linestyle",line_style_model_only,"DisplayName","only background")
       shaded_error_bar((y_pixel_range_ganymede_centred_subimage-y_pixel_center)/R_ganymede_pixel,ganymede_centred_subimage_vertical_box_plot,error_bar_vertical_box_plot,Alpha=shaded_error_alphas,Color=error_bar_color)
       boxing_ax1.XLabel.String = "y-direction [R_{G}]";
       boxing_ax1.YLabel.String = "Brightness  [kR]";
       boxing_ax1.Title.String = observation.name;
%        legend_ax1 = formal_legend(legend("location","southeast"))
fontsize(gca,scale=1.5)     
    end
    plot(boxing_ax1,(y_pixel_range_ganymede_centred_subimage-y_pixel_center)/R_ganymede_pixel,model_image_vertical_box_plot,"color","k","linestyle",line_style_array(counter2 ),"DisplayName",sprintf("n_0 = %.2f [1/cm^3]",n0_value))
    
    
    if counter2 == 1
       boxing_ax2 = formal_axes(subplot(3,3,counter+3));
       hold on
       plot(boxing_ax2,horizontal_range/R_ganymede_pixel,ganymede_centred_subimage_horizontal_box_plot,"DisplayName","STIS observation")
       plot(boxing_ax2,horizontal_range/R_ganymede_pixel,bkg_horizontal_box_plot,"color","k","linestyle",line_style_model_only,"DisplayName","only background")
       shaded_error_bar((x_pixel_range_ganymede_centred_subimage-x_pixel_center)/R_ganymede_pixel,ganymede_centred_subimage_horizontal_box_plot,error_bar_horizontal_box_plot,Alpha=shaded_error_alphas,Color=error_bar_color)
       boxing_ax2.XLabel.String = "x-direction [R_G]";
       boxing_ax2.YLabel.String = "Brightness  [kR]";
%        legend_ax2 = formal_legend(legend("location","southeast"))
        fontsize(gca,scale=1.5)    
    end
    plot(boxing_ax2,(x_pixel_range_ganymede_centred_subimage-x_pixel_center)/R_ganymede_pixel,model_image_horizontal_box_plot,"color","k","linestyle",line_style_array(counter2 ),"DisplayName",sprintf("n_0 = %4.0f [1/cm^3]",n0_value))
    
    if counter2 ==1
        boxing_ax3 = formal_axes(subplot(3,3,counter+6));
        hold on
        plot(boxing_ax3,radial_range/R_ganymede_pixel,ganymede_centred_subimage_radial_plot,"DisplayName","STIS observation")
        plot(boxing_ax3,radial_range/R_ganymede_pixel,bkg_radial_plot,"color","k","linestyle",line_style_model_only,"DisplayName","only background")
        axis([0 2.4 0 1.1*max(bkg_radial_plot) ])
        shaded_error_bar(radial_range/R_ganymede_pixel,ganymede_centred_subimage_radial_plot,error_bar_radial_plot,Alpha=shaded_error_alphas,Color=error_bar_color)
        boxing_ax3.XLabel.String = "Radial distance [R_G]";
        boxing_ax3.YLabel.String = "Brightness  [kR]";
        legend_ax3 = formal_legend(legend("location","southeast"));
        xline(1,':k','HandleVisibility','off')
        fontsize(gca,scale=1.5)    
    end
    plot(boxing_ax3,radial_range/R_ganymede_pixel,model_image_radial_plot,"color","k","linestyle",line_style_array(counter2 ),"DisplayName",sprintf("n_0 = %i [1/cm^3]",n0_value)) 
    
    counter2 = counter2 +1;

end
counter = counter +1;

end
%%
xcorner = [x_pixel_range_ganymede_centred_subimage(1),x_pixel_range_ganymede_centred_subimage(end)];
ycorner  = [y_pixel_range_ganymede_centred_subimage(1),y_pixel_range_ganymede_centred_subimage(end)];

ax = formal_axes(subplot(3,3,3));
hold on
%colormap gray
imagesc(xcorner,ycorner,final_model_image_PSF);
xvertces_vertical_box = [x_pixel_center-globalParameters.diameter_ganymede/2,x_pixel_center-globalParameters.diameter_ganymede/2,x_pixel_center+globalParameters.diameter_ganymede/2,x_pixel_center+globalParameters.diameter_ganymede/2];
yvertces_vertical_box = [y_pixel_range_ganymede_centred_subimage(1),y_pixel_range_ganymede_centred_subimage(end),y_pixel_range_ganymede_centred_subimage(end),y_pixel_range_ganymede_centred_subimage(1)];
fill(xvertces_vertical_box,yvertces_vertical_box,"r",FaceAlpha=0.3)
  axis square;
ax.XLim = [xcorner(1),xcorner(end)];
ax.YLim = [ycorner(1),ycorner(end)];
ax.XLabel.String = "x pixel";
ax.YLabel.String = "y pixel";
fontsize(gca,scale=1.5)    

ax = formal_axes(subplot(3,3,6));
hold on
imagesc(xcorner,ycorner,final_model_image_PSF);
yvertces_vertical_box = [y_pixel_center-globalParameters.diameter_ganymede/2,y_pixel_center-globalParameters.diameter_ganymede/2,y_pixel_center+globalParameters.diameter_ganymede/2,y_pixel_center+globalParameters.diameter_ganymede/2];
xvertces_vertical_box = [x_pixel_range_ganymede_centred_subimage(1),x_pixel_range_ganymede_centred_subimage(end),x_pixel_range_ganymede_centred_subimage(end),x_pixel_range_ganymede_centred_subimage(1)];
fill(xvertces_vertical_box,yvertces_vertical_box,"r",FaceAlpha=0.3)
  axis square;
ax.XLim = [xcorner(1),xcorner(end)];
ax.YLim = [ycorner(1),ycorner(end)];
ax.XLabel.String = "x pixel";
ax.YLabel.String = "y pixel";
fontsize(gca,scale=1.5)    

ax =formal_axes(subplot(3,3,9));
hold on
  cscale                   = [0,7];
imagesc(xcorner,ycorner,final_model_image_PSF,cscale);

theta = linspace(0,2*pi,100);
yvertces_vertical_box = y_pixel_center+globalParameters.diameter_ganymede/2*2.3*cos(theta);
xvertces_vertical_box = x_pixel_center+globalParameters.diameter_ganymede/2*2.3*sin(theta);
%fill(xvertces_vertical_box,yvertces_vertical_box,"r",FaceAlpha=0.3,EdgeColor="non")
  axis square;
ax.XLim = [xcorner(1),xcorner(end)];
ax.YLim = [ycorner(1),ycorner(end)];

ax.XLabel.String = "x pixel";
ax.YLabel.String = "y pixel";

 color_bar_instance = formal_colorbar(colorbar());
    color_bar_instance.Label.String = "Brightness [kR]";
    color_bar_instance.Ticks = 0:2:30;

fontsize(gca,scale=1.5)    
 
%set(gca,'FontSize',14)
%set(findall(gcf,'type','text'),'FontSize',16)
%set(findall(gca,'type','text'),'FontSize',16)



    
    
