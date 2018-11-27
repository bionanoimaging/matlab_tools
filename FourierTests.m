% FourierTests.m : Tests various Fourier transforms for validity and does some speed measurements

a1m = rand(1,100);
a2m = [1 2 3 4 5 6 7 8; 2 3 2 3 2 3 2 3; 4 3 2 4 3 2 4 3; 6 7 6 5 4 5 4 2];
a3m = rand(20,22,26);

a1f = fft(a1m);
a2f = fft(a2m);
a3f = fft(a3m);

a1r = rft(a1m);
a2r = rft(a2m);
a3r = rft(a3m);

GenericDebugCheck('fft','cuda',@castToMatlab,cuda(a1m));
GenericDebugCheck('fft','cuda',@castToMatlab,cuda(a2m));
GenericDebugCheck('fft','cuda',@castToMatlab,cuda(a3m));
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a1m));
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a2m));
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m));

GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a2m),[1 0]);
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a2m),[0 1]);  
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a2m),[1 1]);  
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m),[1 0 0]);
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m),[0 1 0]);
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m),[0 0 1]);
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m),[1 1 0]);
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m),[1 0 1]);
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m),[0 1 1]);
GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m),[1 1 1]);

GenericDebugCheck('ifft','cuda',@castToMatlab,cuda(a1f));
GenericDebugCheck('ifft','cuda',@castToMatlab,cuda(a2f));
GenericDebugCheck('ifft','cuda',@castToMatlab,cuda(a3f));
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a1f));
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a2f));
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a3f));

GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a2r),[1 0]);
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a2r),[0 1]);
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a3r),[1 0 0]);
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a3r),[0 1 0]);
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a3r),[0 0 1]);
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a3r),[1 1 0]);
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a3r),[1 0 1]);
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a3r),[0 1 1]);
GenericDebugCheck('rift','cuda',@castToMatlab,cuda(a3r),[1 1 1]);

% conv1 = real(ifftn(fftn(a3m) .* fftn(a3m)));
% conv2 = real(rift(rft(a3m) .* rft(a3m)));
% dd conv1 conv2 conv1-conv2

%% check some DipImage stuff

a3d = readim('chromo3d');
GenericDebugCheck('ft','cuda',@castToMatlab,cuda(a3d));
GenericDebugCheck('ft1d','cuda',@castToMatlab,cuda(a3d));
GenericDebugCheck('ft2d','cuda',@castToMatlab,cuda(a3d));

a4d = repmat(a3d,[1 1 1 4]);
GenericDebugCheck('ft3d','cuda',@castToMatlab,cuda(a4d));

%% Speed tests
% cuda_shutdown;setCudaSynchronize(1);initCuda();
a2d = readim('orka');
a3d = readim('chromo3d');
a3db = repmat(a3d,[4 4 4]);
a3m = double(a3db);

[TimeInput,TimeCast,res]=GenericDebugCheck('ft3d','cuda',@castToMatlab,cuda(a3d));
[ITimeInput,ITimeCast,res]=GenericDebugCheck('ift3d','cuda',@castToMatlab,cuda(res));
fprintf('\n\n ft3d Cuda %g, DipImage %g, ift3d Cuda %g, DipImage %g, \n',TimeInput,TimeCast,ITimeInput,ITimeCast)

[TimeInput,TimeCast,res]=GenericDebugCheck('ft2d','cuda',@castToMatlab,cuda(a2d));
[ITimeInput,ITimeCast]=GenericDebugCheck('ift2d','cuda',@castToMatlab,cuda(res));
fprintf('2D Data ft2d Cuda %g, DipImage %g, ift2d Cuda %g, DipImage %g, \n',TimeInput,TimeCast,ITimeInput,ITimeCast)

[TimeInput,TimeCast,res]=GenericDebugCheck('rft3d','cuda',@castToMatlab,cuda(a3d));
[ITimeInput,ITimeCast]=GenericDebugCheck('rift3d','cuda',@castToMatlab,cuda(res));
fprintf('3D Data rft3d Cuda %g, DipImage %g, rift3d Cuda %g, DipImage %g, \n',TimeInput,TimeCast,ITimeInput,ITimeCast)

[TimeInput,TimeCast,res]=GenericDebugCheck('fftn','cuda',@castToMatlab,cuda(a3m));
[ITimeInput,ITimeCast]=GenericDebugCheck('ifftn','cuda',@castToMatlab,cuda(res));
fprintf('Large 3D Data fftn Cuda %g, Matlab %g, ifftn Cuda %g, Matlab %g, \n',TimeInput,TimeCast,ITimeInput,ITimeCast)

% profile clear; profile on -history
[TimeInput,TimeCast,res]=GenericDebugCheck('rft','cuda',@castToMatlab,cuda(a3m));
[ITimeInput,ITimeCast]=GenericDebugCheck('rift','cuda',@castToMatlab,cuda(res));
fprintf('Large 3D Data rft Cuda %g, Mex %g, rift Cuda %g, Mex %g, \n',TimeInput,TimeCast,ITimeInput,ITimeCast)
% profile off; profview
