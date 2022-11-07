function [fitresult, gof] = Fourier(X, I)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( X, I );

% Set up fittype and options.
ft = fittype( 'fourier5' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% AMPS=[fitresult.b5 fitresult.b4 fitresult.b3 fitresult.b2 fitresult.b1 fitresult.a0 fitresult.a1 fitresult.a2 fitresult.a3 fitresult.a4 fitresult.a5]
% bar([-5 -4 -3 -2 -1 0 1 2 3 4 5], AMPS)

