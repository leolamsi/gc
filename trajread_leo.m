function [m] = trajread_leo(filename, T)%, timestep)
    % developed for glassy crystal project by Leo
    % filename: trajectory file name, T: temperature (eg 0.7), timestep:
    % timestep used in lammps simulation
    % reads lammps dump: (dump 2 all custom ${dumpframe} ${out_traj_path} id x y z)
    % read trajectory, unwrap periodic boundary, calculate msd,
    % relaxation time tau and diffusion coefficient D
    
    % Open the file
    fid = fopen(filename, 'r');
    % Check if the file opened successfully
    if fid == -1
        error('Cannot open the file.');
    end
    
    % Initialize variables
    atomData = [];
    timesteps = [];
    numAtoms = 0;
    boxBounds = [];
    m.nframe = 0;
    m.T = T;
    
    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        
        if contains(line, 'ITEM: TIMESTEP')
            m.nframe = m.nframe + 1;
            % Read the timestep
            timesteps = [timesteps; str2double(fgetl(fid))];
        elseif contains(line, 'ITEM: NUMBER OF ATOMS')
            % Read the number of atoms
            numAtoms = str2double(fgetl(fid));
        elseif contains(line, 'ITEM: BOX BOUNDS')
            % Read the box bounds
            boxBounds = zeros(3, 2);
            for i = 1:3
                boxBounds(i, :) = str2double(strsplit(fgetl(fid)));
            end
        elseif contains(line, 'ITEM: ATOMS id x y z')
            % Read the atom data
            tempData = zeros(numAtoms, 4);
            for i = 1:numAtoms
                tempData(i, :) = str2double(strsplit(strtrim(fgetl(fid))));
            end
            atomData(m.nframe, :, :) = tempData(:, 2:end);
        end
    end
    
    % Close the file
    fclose(fid);
    
    % Assign output structure
    m.d = 3;
    m.r = atomData;
    m.dt = timesteps(2)*0.001;% * timestep; % convert from timesteps to LJ time
    m.n = numAtoms;
    m.L = boxBounds(1, 2);
    m.rmin = zeros(m.d, 1);
    m.rmax = ones(m.d, 1) * m.L;


    % unwrap
    L = m.rmax - m.rmin;
    for fr = 2:m.nframe
        for k = 1:m.d
            dr = m.r(fr,:,k) - m.r(fr-1,:,k);
            m.r(fr,:,k) = m.r(fr,:,k) - L(k) * round(dr/L(k));
        end
    end
    
    m = msdcal_leo(m, T, filename(1:end-5));
    fscal(m,T, filename(1:end-5));
end