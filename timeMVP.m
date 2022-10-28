% timeMVP.m

function varargout = time (varargin)

    [vararout{1:nargout}] = feval(varargin{:});

end

% handles KVFs for all values of N, tolerances 1e-6 and 1e-8
function r = build()
  addpaths();

  Ns = [4, 50, 100, 500, 1000];
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
    c.X = [X; Y; Z];

    %  create kernel based on mesh
    G = 30e3;
    nu = 0.25; % this is fine
    % theta is zero
    % eps12 can be 1; others zero

    disp('mesh done, making kernel')

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

    disp('kernel done')

    % general things
    c.command = 'compress';
    c.greens_fn = 'time';
    c.err_method = 'mrem-fro';
    %c.K = real(kernel(:)');
    %size(c.K)
    c.allow_overwrite = 1;

    % tolerance 1e-6;
    c.tol = 1e-6;
    %c.write_hmat_filename = strcat('./tmp/t3d_n',nstring,'_e6');
    c.write_hmat_filename = './tmp/timingHM';
    %c.write_hd_filename = './tmp/timingHM-hd';
    c.kvf = [c.write_hmat_filename, '.kvf'];
    c.kvf
    kvf('Write', c.kvf, c, 4);
    cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
    disp(cmd)

    % tolerance 1e-8;
    c.tol = 1e-8;
    c.write_hmat_filename = ['./tmp/t3d_n' nstring '_e8'];
    c.write_hd_filename = [c.write_hmat_filename '-hd'];
    c.kvf = [c.write_hmat_filename '.kvf'];
    kvf('Write', c.kvf, c, 32);
    cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
    disp(cmd)
    %r.kvfs_2d_8(length(r.kvfs_3d_8)+1) = c.write_hmat_filename;

  end

  b.Ns = Ns;

end

function measure(b)

  Ns = b.Ns;
  t_3d_6 = [];
  t_3d_8 = [];
  t_dense = [];

  Ns = [10, 25, 50, 100, 500, 1000];
  for i=1:length(Ns)
    N = Ns(i);
    v = ones(N*N*N)

    tees = [];
    for j=1:3
      hme6 = hmmvp('load', './tmp/timing_3d_n' + nstring + '_e6');
      now = tic();
      x = hmmvp('mvp', hme6, v);
      t = toc(now)
      tees(j) = t;
    end
    t_3d_6(i) = mean(tees);

    tees = [];
    for j=1:3
      hme8 = hmmvp('load', './tmp/timing_3d_n' + nstring + '_e8');
      now = tic();
      x = hmmvp('mvp', hme8, v);
      t = toc(now)
      tees(i) = t;
    end
    t_3d_8(i) = mean(tees);

    tees = [];
    for j = 1:3
      dense = hmmvp('extract', hme8, 1:1:N*N*N, 1:1:N*N*N);
      now = tic();
      x = dense*v;
      t = toc(now)
      tees(i) = t;
    end
    t_dense(i) = mean(tees);

  end

  T = [t_dense; t_3d_6; t_3d_8]

end

% XXX --- Private --- XXX
function addpaths()
  addpath('include/hmmvp/matlab')
  addpath('ODEsolving')
  addpath('include/')
  addpath('./tmp/')
end
