function plot_comparison_current
    clc;
    clear all;

    filename1 = 'OutputData - Sim 2 - No Interactions\current_vs_time.dat';
    filename2 = 'OutputData - Sim 1 - Interactions\current_vs_time.dat';
    filename3 = 'OutputData - Sim 2 - No Interactions\coverage_with_analytical_solution.dat';

    data1 = load(filename1);
    stop_array1 = length(data1(:,1));
    for idx1 = 2:length(data1(:,1))
        for idx2 = 4:length(data1(1,:))
            if data1(idx1,idx2) < 1.0e-12 && idx1 < stop_array1
                stop_array1 = idx1;
            end
        end
    end

    data2 = load(filename2);
    stop_array2 = length(data2(:,1));
    for idx1 = 2:length(data2(:,1))
        for idx2 = 4:length(data1(2,:))
            if data1(idx1,idx2) < 1.0e-12 && idx1 < stop_array2
                stop_array2 = idx1;
            end
        end
    end

    time1 = data1(1:stop_array1,2);
    time2 = data2(1:stop_array2,2);

    analytical_time_end = time1(end);
    delta_analytical_time = analytical_time_end/length(time1);
    analytical_time = delta_analytical_time:delta_analytical_time:(analytical_time_end);

    data3 = load(filename3);
    analytical_site_coverage = data3(1:stop_array1,3);
    r_filename = 'reaction_rates.dat';
    r_data = load(r_filename);
    ra = r_data(1,2);%1.4e-7; (site*sec)^(-1)
    rd = r_data(2,2); %2.8e-7; (site*sec)^(-1)

    analytical_sites_to_desorb = analytical_site_coverage; %.*(25*25);
    analytical_sites_desorbed = rd.*analytical_sites_to_desorb;
    analytical_electrons = analytical_sites_desorbed*4;
    analytical_charge_density = analytical_electrons.*(1.6e-19/1.039e-16); %C/m^2

%     analytical_current_density = ((1.0-analytical_sim).*((4*1.6e-19)/1.039e-16))./time1;
    analytical_current_density = zeros(size(analytical_charge_density));
    for idx = 1:length(analytical_charge_density)
        analytical_current_density(idx) = analytical_charge_density(idx)*analytical_time(idx);
    end

%     analytical_sim = data(:,3);
    kmc_time1 = data1(:,2:2:20);
    kmc_time2 = data2(:,2:2:20);

    kmc_sim1 = data1(:,3:2:21);
    kmc_sim2 = data2(:,3:2:21);

    kmc_avg_time_1 = zeros(size(time1));
    kmc_avg_time_2 = zeros(size(time2));

    kmc_avg_current_1 = zeros(size(time1));
    kmc_avg_current_2 = zeros(size(time2));

    for idx = 1:length(time1)
        kmc_avg_time_1(idx) = mean(kmc_time1(idx,:));
        kmc_avg_current_1(idx) = mean(kmc_sim1(idx,:)).*(kmc_avg_time_1(idx)*((4*1.6e-19)/1.039e-16)); %kmc_avg_time_1(idx);
    end

    for idx = 1:length(time2)
        kmc_avg_time_2(idx) = mean(kmc_time2(idx,:));
        kmc_avg_current_2(idx) = mean(kmc_sim2(idx,:)).*(kmc_avg_time_2(idx)*((4*1.6e-19)/1.039e-16));  %kmc_avg_time_2(idx);
    end

    % =====================================================================
    % Plotting variables
    tick_label_size = 20;
    axis_label_size = 24;
    title_label_size = 20;
    plot_line_width = 3;
    axis_line_width = 3;
    marker_size = 2;
    font_weight = 'bold';
    % =====================================================================

    figure(3)
    hold on
%     plot(time, theta_e_plot, '--k','LineWidth', plot_line_width)
%     plot(time,analytical_sim,'-b','LineWidth', plot_line_width)
    plot(kmc_avg_time_1,kmc_avg_current_1.*1.0e4,'-b','LineWidth', plot_line_width)
    plot(kmc_avg_time_2,kmc_avg_current_2.*1.0e4,'-r','LineWidth', plot_line_width)
    plot(analytical_time,analytical_current_density.*1.0e4,'-k','LineWidth', plot_line_width)
%     plot(kmc_time,kmc_sim,'g^','MarkerSize',marker_size)

    axis square
    box on
    xlabel('Time (s)')
    ylabel('Average current density (A/cm^2)')

    set(gca,'FontSize',tick_label_size,'FontWeight',font_weight, ...
    'LineWidth',axis_line_width,'XMinorTick','on','YMinorTick','on','yscale','log')

    legend('KMC model - no interactions', 'KMC model - site blocking', 'Analytical model','location','southeast')
    legend boxoff
    hold off


end

