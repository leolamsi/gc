function [m] = msdcal_leo(m, T, filename)
    % msd
    frame = 1:m.nframe; % all frames
    framestart = frame(1); % ref. time frame
    frameend = min(m.nframe, frame(end));
    sampling_interval = 1000;
    fr_increment = 1.4;
    
    fr_size=1;
    dfrcnt=1;
    dim = 1:m.d;
    r2mean = NaN(m.nframe);
    dt = NaN(m.nframe);
    while fr_size < frameend-framestart 
        dt(dfrcnt) = m.dt*fr_size;
        framesample = frameend-fr_size : -sampling_interval : framestart;
        r2_sum = 0;
        for i = 1:length(framesample)
            fr = framesample(i);
    	    dr = squeeze( m.r(fr+fr_size,:,:) - m.r(fr,:,:) ); % disp vector of all atoms
            dr2 = dot(dr(:,dim), dr(:,dim), 2);
            r2_sum = r2_sum + double(mean( dr2 ))/length(framesample); % average over all atoms, in dim
        end
        r2mean(dfrcnt) = r2_sum; % average over all initial frames
        dfrcnt = dfrcnt+1;
        %fr_increment = 1;
        fr_size = round(fr_size*fr_increment+1); % factor determines density of points in t axis
    end
    r2mean = r2mean(~isnan(r2mean));
    dt = dt(~isnan(dt));

    data = [dt r2mean]; % time and msd
    fname = filename + ".msd";
    save('-ascii', fname, 'data');
    
    D_from_msd(T, filename);
end