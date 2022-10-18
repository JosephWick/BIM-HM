% timeMVP.m

function varargout = time (varargin)

    [vararout{1:nargout}] = feval(varargin{:});

end

% handles KVFs for all values of N, tolerances 1e-6 and 1e-8
function r = build()
  addpaths();

  Ns = [10, 25, 50, 100, 500, 1000];
  probDim = 200e3; %200km
  transition = 40e3;

  r.kvfs_2d_6 = [];
  r.kvfs_2d_8 = [];
  r.kvfs_3d_6 = [];
  r.kvfs_3d_8 = [];
  r.N = Ns;

  disp('Run the following commands in a shell to build all HMs')

  for i=1:length(Ns)
    N = Ns(i);
    nstring = string(N);

    eps = 1e-12;
    nc = (-N/2:N/2);

    Xhat = linspace(0, probDim, N);
    Yhat = tan(nc*pi/(2.5*max(nc)))*32e3;
    Zhat = transition+tan((0:N)'*pi/(2.2*(N+eps)))*transition;

    L = abs(Xhat(2) - Xhat(1));

    % --- 3D Mesh ---
    % mesh should not be unifornm
    [Z, Y, X] = ndgrid(Zhat, Yhat, Xhat);
    X = X(:)';
    Y = Y(:)';
    Z = Z(:)';
    c.X = [X(:)', Y(:)', Z(:)'];

    %  create kernel based on mesh
    G = 30e3;
    nu = 0.25; % this is fine
    % theta is zero
    % eps12 can be 1; others zero

    kernel = zeros(N^3, N^3);
    for i = 1:N^3
      for j = 1:N^3
        kernel(i,j) = computeStressVerticalShearZone_s12(...
        X(i), Y(i), Z(i), ...
        X(j), Y(j), Z(j), L, L, L, 0, ...
        0, 1, 0, 0, 0, 0, ...
        G, nu);

      end
    end

    % general things
    c.greens_fn = 'shear1212';
    c.command = 'compress';
    c.err_method = 'mrem-fro';
    c.K = kernel;

    % tolerance 1e-6;
    c.tol = 1e-6;
    c.write_hmat_filename = './tmp/timing_3d_n' + nstring + '_e6';
    c.write_hd_filename = c.write_hmat_filename + '-hd';
    c.kvf('Write', c.kvf, c, 32);
    cmd = '    include/hmmvp/bin/hmmvpbuild_omp' + c.kvf;
    disp(cmd)
    r.kvfs_3d_6(length(r.kvfs_3d_6)+1) = c.write_hmat_filename;

    % tolerance 1e-8;
    c.tol = 1e-8;
    c.write_hmat_filename = './tmp/timing_3d_n' + nstring + '_e8';
    c.write_hd_filename = c.write_hmat_filename + '-hd';
    c.kvf = c.write_hmat_filename + '.kvf'
    kvf('Write', c.kvf, c, 32);
    cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
    disp(cmd)
    r.kvfs_2d_8(length(r.kvfs_3d_8)+1) = c.write_hmat_filename;

  end

end

function measure(b)

  Ns = b.Ns;
  t_2d_6 = [];
  t_2d_8 = [];
  t_3d_6 = [];
  t_3d_8 = [];



end

% XXX --- Private --- XXX
function addpaths()
  addpath('include/hmmvp/matlab')
  addpath('ODEsolving')
  addpath('include/')
end
