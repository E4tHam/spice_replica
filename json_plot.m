
clear;

% read json into struct
filename = 'out.json';
tran = jsondecode(fileread(filename));
fclose('all');

stop_time = ceil(tran.stop_time / tran.time_step) * tran.time_step;
x_ticks = (0):(tran.time_step):(stop_time);

for i = 1:length(tran.NODES)
    for j = 1:length(tran.PLOTNV)
        if tran.NODES(i).name == tran.PLOTNV(j)
            figure();
            plot(x_ticks, tran.NODES(i).voltages);
            title(append('Node ', string(tran.NODES(i).name), ' voltage'));
            xlabel('time (s)');
            ylabel('Voltage');
            xlim([0 tran.stop_time]);
            ylim([-inf inf]);
            if all(tran.NODES(i).voltages >= 0)
                ylim([0 inf]);
            end
            if all(tran.NODES(i).voltages <= 0)
                ylim([-inf 0]);
            end
        end
    end
end
