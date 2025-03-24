function D_from_msd(T, filename)
    % for 3 dimensions
    ndim = 3;

    cutoff_grad = 0.96; % for D

    g0D = 0.1; % diffusive regime when msd > 1


    % Open the file
    fid = fopen(filename+".msd", 'r');
    
    % Read the data from the file using textscan
    data = textscan(fid, '%f %f', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    
    % Close the file
    fclose(fid);
    dt = data{1};
    r2mean = data{2};

    % find ralexation time tau

    g0tau = 0.05; % msd when still experiencing cage effect
    tau = 0;
    ind = find( r2mean<g0tau, 1, 'last' );%find the last position r2mean<opt.g0tau
    if ~isempty(ind) && ind<size(r2mean,2)
      tau = interp1( r2mean([ind ind+1]), dt([ind ind+1]), g0tau );
    end
    
    % find diffusion coefficient D

    D = 0;
    dto = dt(r2mean>g0D); % points outside cage (diffusive regime)
    r2o = r2mean(r2mean>g0D); % points outside cage
    if size(dto, 2) > 0
        fin=1;
        grad = zeros(length(1:length(r2o)-1));
        for i = 1:length(r2o)-1
            dxlog = log(dto(i+1)/dto(i));
            dylog = log(r2o(i+1)/r2o(i));
            grad(i) = dylog/dxlog;
        end  
        
        for i = 1:length(grad)
            fin = i;
            if grad(i) > cutoff_grad            %cutoff gradient for selecting points in diffusive regime
                break     
            end
        end
        if length(r2o)-fin == 1
            fprintf('Diffusive regime not reached with gradient %d\n', grad(fin));
        else
            fprintf('%d points used, gradient %d\n', length(r2o)-fin, grad(fin));
        end
        
        r2oD = r2o(fin:end);
        dtoD = dto(fin:end);
        D = mean(r2oD./dtoD)/(2*ndim);
        %msd.cut = [dtoD(1), r2oD(1)];
    end
    
    data = [T tau D];
    fname = filename + ".tauD";
    save('-ascii', fname, 'data');
end