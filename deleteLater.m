for i = 1:length(logger)
    for  j = 1:length(logger(i).flow.t)
        qlog = logger(i).flow.states.x(:,j);

        disp.update(qlog);

        pause(0.2);
    end
end
