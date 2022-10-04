% timeMVP.m

function varargout = time (varargin)

    [vararout{1:nargout}] = feval(varargin{:});

end

% handles KVFs for all values of N, tolerances 1e-6 and 1e-8
function r = build()
  addpaths();

  Ns = [10, 25, 50, 100, 500, 1000];
  probDim = 200e3; %200km

  c.greens_fn = 'shear1212';

  r.kvfs_2d_6 = [];
  r.kvfs_2d_8 = [];
  r.kvfs_3d_6 = [];
  r.kvfs_3d_8 = [];
  r.N = Ns;

  disp('Run the following commands in a shell to build all HMs')

  for i=1:length(Ns)
    N = Ns(i);
    nstring = string(N);

    Xhat = linspace(0, probDim, N);
    Yhat = linspace(0, probDim, N);
    Zhat = linspace(0, probDim, N);

    % --- 2D mesh ---
    [Y, X] = ndgrid(Yhat, Xhat);
    Z = zeros(size(X(:)'));

    % --- 2D kvf ---
    c.X = [X(:)', Y(:)', Z];

    % tolerance 1e-6;
    c.tol = 1e-6;
    c.write_hmat_filename = strcat('./tmp/timing_2d_n', nstring, '_e6');
    c.write_hd_filename = strcat(c.write_hmat_filename, '-hd');
    c.kvf = strcat(c.write_hmat_filename, '.kvf');
    kvf('Write', c.kvf, c, 32);
    cmd = '    include/hmmvp/bin/hmmvpbuild_omp ' + c.kvf;
    disp(cmd)
    r.kvfs_2d_6(length(r.kvfs_2d_6)+1) = c.write_hmat_filename;

    % tolerance 1e-6;
    c.tol = 1e-8;
    c.write_hmat_filename = './tmp/timing_2d_n' + nstring + '_e8';
    c.write_hd_filename = c.write_hmat_filename + '-hd';
    c.kvf = c.write_hmat_filename + '.kvf';
    kvf('Write', c.kvf, c, 32);
    cmd = '    include/hmmvp/bin/hmmvpbuild_omp ' + c.kvf;
    disp(cmd)
    r.kvfs_2d_8(length(r.kvfs_2d_8)+1) = c.write_hmat_filename;

    % --- 3D Mesh ---
    [Z, Y, X] = ndgrid(Zhat, Yhat, Xhat);
    c.X = [X(:)', Y(:)', Z(:)'];

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
