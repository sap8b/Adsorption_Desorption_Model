function plot_theta
    clc;
    clear all;

    filename = 'coverage_with_analytical_solution.dat';
    data = load(filename);

    time = data(:,2);
    analytical_sim = data(:,3);
    kmc_time = data(:,4:2:22);
    kmc_sim = data(:,5:2:23);

    kmc_avg_time = zeros(size(time));
    kmc_avg_coverage = zeros(size(time));

    for idx = 1:length(time)
        kmc_avg_time(idx) = mean(kmc_time(idx,:));
        kmc_avg_coverage(idx) = mean(kmc_sim(idx,:));
    end

    ra = 1.4e-7;
    rb = 2.8e-7;
    theta_e = ra/(ra+rb);
    theta_e_plot = zeros(size(time));
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
    plot(time, theta_e_plot, '--k','LineWidth', plot_line_width)
    plot(time,analytical_sim,'-b','LineWidth', plot_line_width)
    plot(kmc_avg_time,kmc_avg_coverage,'-r','LineWidth', plot_line_width)
%     plot(kmc_time,kmc_sim,'g^','MarkerSize',marker_size)

    axis square
    box on
    xlabel('Time (s)')
    ylabel('Surface coverage')

    set(gca,'FontSize',tick_label_size,'FontWeight',font_weight, ...
    'LineWidth',axis_line_width,'XMinorTick','on','YMinorTick','on')

    legend('Equilibrium coverage', 'Analytical solution','KMC model average', 'KMC model run results','location','southeast')
    legend boxoff
    hold off


end