function plot_coverage_sites
    clc;
    clear all;
    filename = 'configuration.dat';

    raw_site_data = load(filename);
    num_atoms = length(raw_site_data(:,4));
    marker_face_colors = [1:num_atoms; 1:num_atoms; 1:num_atoms]';
    sz = zeros(1,num_atoms);

    al_color = [0.5, 0.5, 0.7];
    cr_color = [0.1, 0.1, 0.6];
    o2_color = [0.8, 0.1, 0.5];
    fe_color = [0.3, 0.3, 0.3];
    ag_color = [1.0, 1.0, 1.0];
    cl_color = [0.9, 0.1, 0.1];
    H2O_color = [0.0, 0.0, 0.8];
    OH_color = [0.0, 0.0, 1.0];
    H_Color = [0.0, 0.0, 0.6];
    au_color = [0.9 0.9 0.0];

    for count = 1:num_atoms
        switch raw_site_data(count,4)
            case 0
                marker_face_colors(count,:) = [1 1 28];
                sz(1,count) = 10;
            case 1
                marker_face_colors(count,:) = o2_color;
                sz(1,count) = 50;
             case 2
                marker_face_colors(count,:) = cr_color;
                sz(1,count) = 70;
            case 3
                marker_face_colors(count,:) = al_color;
                sz(1,count) = 50;
            case 4
                marker_face_colors(count,:) = fe_color;
                sz(1,count) = 50;
            case 5
                marker_face_colors(count,:) = ag_color;
                sz(1,count) = 40;
            case 6
                marker_face_colors(count,:) = cl_color;
                sz(1,count) = 90;
            case 7
                marker_face_colors(count,:) = H2O_color;
                sz(1,count) = 1;
            case 8
                marker_face_colors(count,:) = OH_color;
                sz(1,count) = 10;
            case 9
                marker_face_colors(count,:) = au_color;
                sz(1,count) = 70;
        end
    end
    figure(2)
    hold on
%     scatter3(raw_atom_data(:,1),raw_atom_data(:,2),raw_atom_data(:,3), ...
%         sz, ...
%         'MarkerEdgeColor',[0 0 0], ...
%         'MarkerFaceColor',marker_face_colors(:,:), ...
%         'LineWidth',1.5)
    h = scatter3(raw_site_data(:,1),raw_site_data(:,2),raw_site_data(:,3), ...
        sz, ...
        marker_face_colors, ...
        'filled');
    h.MarkerEdgeColor = [0 0 0];
%     axis square
    box on
    xlabel('x position (nm)')
    ylabel('y position (nm)')
    zlabel('z position (nm)')
    %zlim([0 0.01*max(raw_site_data(:,1))])

    az = 90;
    el = 90;
    view(az, el);

    hold off
end
