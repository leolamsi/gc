function [fs] = fscal(m, T, filename, varargin)
% modified from fpcfcal.m
% Calculate sisf (self-intermediate scattering function) and
%			overlap_q (local particle persistence) and
%           fpcf (four-point correlation function)

% Output:
%   fpcf.t, fpcf.Fs(ik,idt), fpcf.overlap_q
% Arguments:
%   k        : reduced wave vector for mobility
%   q        : reduced wave vector for Fourier
%   interval : sampling interval
%   frinc    : fractional increament of dt

opt = getopt(struct('lambda', [2], ...
                    'q', [0], ...
                    'interval', 10, ...
                    'frinc', 1.4), varargin{:});

%k = 2 * pi / m.L * opt.k(:); nk = length(k);
k = 2 * pi / opt.lambda(:); nk = length(k);
q = 2 * pi / m.L * opt.q(:); nq = length(q);

dfx = 1; idt = 1;
while dfx < m.nframe
    df(idt) = dfx;
    fpcf.t(idt) = m.dt * dfx;
    idt = idt + 1;
    dfx = round(dfx * opt.frinc + 1);
end
ndt = length(df);
% Fs
for idt = 1:ndt
    frsamp = m.nframe-df(idt) : -opt.interval : 1;
    dr = m.r(frsamp+df(idt),:,:) - m.r(frsamp,:,:);

    %overlap_q
    oq = dr(:)';
    oq = oq == 0;
    fpcf.overlap_q(:,idt) = mean(oq,2);
    %overlap_q

    fpcf.Fs(:,idt) = mean(cos(k*dr(:)'), 2);
end

data = [fpcf.t', fpcf.Fs'];
fname = filename + ".fs";
save('-ascii', fname, 'data');

end
