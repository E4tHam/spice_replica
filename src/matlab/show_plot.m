
function show_plot(in_data, in_title, in_xlabel, in_ylabel, in_xtick, in_xlim)

    x_ticks = (0):(length(in_data)-1);
    x_ticks = x_ticks * in_xtick;

    figure();
    plot(x_ticks, in_data);
    title(in_title);
    xlabel(in_xlabel);
    ylabel(in_ylabel);
    xlim([0 in_xlim]);
    ylim([-inf inf]);
    if all(in_data >= 0)
        ylim([0 inf]);
    end
    if all(in_data <= 0)
        ylim([-inf 0]);
    end

end
