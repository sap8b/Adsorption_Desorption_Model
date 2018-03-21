function plot_comparison_theta
    clc;
    clear all;

    filename1 = 'OutputData - Sim 2 - No Interactions\coverage_with_analytical_solution.dat';
    filename2 = 'OutputData - Sim 1 - Interactions\coverage_with_analytical_solution.dat';
    data1 = load(filename1);
    data2 = load(filename2);

    time1 = data1(:,2);
    time2 = data2(:,2);
    analytical_sim = data1(:,3);

%     analytical_sim = data(:,3);
    kmc_time1 = data1(:,4:2:22);
    kmc_time2 = data2(:,4:2:22);

    kmc_sim1 = data1(:,5:2:23);
    kmc_sim2 = data2(:,5:2:23);

    kmc_avg_time_1 = zeros(size(time1));
    kmc_avg_time_2 = zeros(size(time2));

    kmc_avg_coverage_1 = zeros(size(time1));
    kmc_avg_coverage_2 = zeros(size(time2));

    for idx = 1:length(time1)
        kmc_avg_time_1(idx) = mean(kmc_time1(idx,:));
        kmc_avg_coverage_1(idx) = mean(kmc_sim1(idx,:));
    end

    for idx = 1:length(time2)
        kmc_avg_time_2(idx) = mean(kmc_time2(idx,:));
        kmc_avg_coverage_2(idx) = mean(kmc_sim2(idx,:));
    end

    r_filename = 'reaction_rates.dat';
    r_data = load(r_filename);
    ra = r_data(1,2);%1.4e-7;
    rb = r_data(2,2); %2.8e-7;
    theta_e = ra/(ra+rb);
    theta_e_plot = zeros(size(time1));
    theta_e_plot(1:end) = theta_e;

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

    figure(1)
    hold on
    plot(time1, theta_e_plot, '--k','LineWidth', plot_line_width)
    plot(time1,analytical_sim,'-b','LineWidth', plot_line_width)
    plot(kmc_avg_time_1,kmc_avg_coverage_1,'-r','LineWidth', plot_line_width)
    plot(kmc_avg_time_2,kmc_avg_coverage_2,'-k','LineWidth', plot_line_width)
% %     plot(kmc_time,kmc_sim,'g^','MarkerSize',marker_size)

    axis square
    box on
    xlabel('Time (s)')
    ylabel('Surface coverage')

    set(gca,'FontSize',tick_label_size,'FontWeight',font_weight, ...
    'LineWidth',axis_line_width,'XMinorTick','on','YMinorTick','on')

    legend('Equilibrium coverage', 'Analytical solution coverage','KMC model predicted coverage - No interactions', 'KMC model predicted coverage - Site blocking','location','southeast')
    legend boxoff
    hold off


end