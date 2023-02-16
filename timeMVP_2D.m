% timeMVP.m

function varargout = tMVP (varargin)

    [varargout{1:nargout}] = feval(varargin{:});

end

% handles KVFs for all values of N, tolerances 1e-6 and 1e-8
function r = build()
  addpaths();

  Ns = [30, 90, 164, 250, 269];
  % total elems in mesh: [900, 8100, 27000, 62500, 75,600]
  probDim = 200e3; %200km
  transition = 40e3;

  disp('Run the following commands in a shell to build all HMs')

  i = 5;
  N = Ns(i);
  r.N = N;

  eps = 1e-12;
  nc = (-N/2:N/2);

  Yhat = tan(nc*pi/(2.5*max(nc)))*32e3;
  Zhat = transition+tan((0:N)'*pi/(2.2*(N+eps)))*transition;
  Xhat = zeros(size(Yhat));

  L = zeros(1,N);
  W = zeros(1,N);

  Xchat = zeros(1,N);
  Ychat = zeros(1,N);
  Zchat = zeros(1,N);
  for idx=(1:length(Zhat)-1)
    L(idx) = abs(Yhat(idx) - Yhat(idx+1));
    W(idx) = abs(Zhat(idx) - Zhat(idx+1));

    Ychat(idx) = Yhat(idx) - abs(Yhat(idx) - Yhat(idx+1))/2;
    Zchat(idx) = Zhat(idx) - abs(Zhat(idx) - Zhat(idx+1))/2;
  end
  L(end) = L(1);
  W(end) = abs(Zhat(end-1) - Zhat(end));
  Yhat(end)=[]; Zhat(end)=[]; Xhat(end)=[];

  length(Xhat);
  length(Yhat);
  length(Zhat);

  % --- Package for KVFs ---
  [shearY shearZ] = ndgrid(Yhat, Zhat);
  [shearY_c shearZ_c] = ndgrid(Ychat, Zchat);

  shearX = zeros(size(shearY));
  shearX_c = zeros(size(shearY_c));

  c.L = L;
  c.W = W;

  c.X = [shearX_c(:)'; shearY(:)'; shearZ(:)'];
  c.Y = [Xhat; Ychat; Zchat];

  c.Z = [shearX(:)'; shearY(:)'; shearZ(:)'];

  c.Ny = N;
  c.Nz = N;

  % create kernel based on mesh
  c.G = 30e3;
  c.nu = 0.25; % this is fine
  % theta is zero
  % eps12 can be 1; others zero

  % general things
  c.command = 'compress';
  c.greens_fn = 'shear1212';
  c.err_method = 'mrem-fro';
  c.allow_overwrite = 1;

  %c.Bfro = 1e-8;

  % tolerance 1e-6;
  c.tol = 1e-6;
  c.write_hmat_filename = './tmp/t2d_e6';
  c.write_hd_filename = './tmp/t2d_e6-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.kvfe6 = c.write_hmat_filename;

  % tolerance 1e-8;
  c.tol = 1e-8;
  c.write_hmat_filename = './tmp/t2d_e8';
  c.write_hd_filename = './tmp/t2d_e8-hd';
  c.kvf = [c.write_hmat_filename '.kvf'];
  kvf('Write', c.kvf, c, 4);
  cmd = ['    include/hmmvp/bin/hmmvpbuild_omp ' c.kvf];
  disp(cmd)
  r.kvfe8 = c.write_hmat_filename;

end

function measure(b)
  N = b.N;
  disp(N)

  t_3d_6 = -1;
  t_3d_8 = -1;
  t_dense = -1;

  rng('default');
  v = zeros(N,N);
  ia = floor(N/4);
  ib = floor(3*N/4);
  v(ia:ib,ia:ib) = 1.0;
  v = v(:);

  numIter = 500;


  disp('timing e6...')
  t6 = [];
  for j=1:5
    hme6 = hmmvp('init', b.kvfe6);
    now = tic();
    x = hmmvp('mvp', hme6, v);
    t = toc(now);
    t6(j) = t;
  end
  t_3d_6 = mean(t6);

  disp('timing e8...')
  t8 = [];
  for j=1:5
    hme8 = hmmvp('init', b.kvfe8);
    now = tic();
    x = hmmvp('mvp', hme8, v);
    t = toc(now);
    t8(j) = t;
  end
  t_3d_8 = mean(t8);

  disp('timing dense...')
  td = [];
  for j = 1:5
    dense = hmmvp('extract', hme8, 1:1:N*N, 1:1:N*N);
    now = tic();
    x = dense*v;
    t = toc(now);
    td(j) = t;
  end
  t_dense = mean(td);

disp('RESULTS:')
disp('dense: ')
t_dense
disp('e6')
t_3d_6
disp('e8')
t_3d_8

end

% XXX --- Private --- XXX
function addpaths()
  addpath('include/hmmvp/matlab')
  addpath('ODEsolving')
  addpath('include/')
  addpath('./tmp/')
end
