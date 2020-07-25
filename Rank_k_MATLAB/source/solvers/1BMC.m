%% 1BMC.m
% 
% Edits by: Mel Beckerleg 2020
% This file was derived from :

% code for a simple demo of the one bit matrix 
% completion software package that supplements the paper 
% http://arxiv.org/abs/1209.3672. 
%
% Most recent change - 5/16/2014
%
% Copyright 2013, M. Davenport, Y. Plan, E. van den Berg, M. Wootters
%
% This file is part of 1BMC Toolbox version 1.2.
%
%    The 1BMC Toolbox is free software: you can redistribute it and/or 
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.
%
%    The 1BMC Toolbox is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with the 1BMC Toolbox. If not, see <http://www.gnu.org/licenses/>.


function [Mhat] = OneBitMC(m,n,r,B,W_omega,opts);
	
	epsilon = opts.epsilon;
	f = @(x) x*(2 e
	f = 
	idx = 

	%% Set up optimization problem
	options = struct();
	options.iterations = 10000; 
	options.stepMax    = 10000;
	options.stepMin    = 1e-4;
	options.optTol     = 1e-3;
	options.stepMax    = 1e9;
		
	funObj  = @(x) logObjectiveGeneral(x,y,idx,f,fprime);

	%% Define alpha to be the correct maximum using an oracle
	alpha   = 1;
	radius  = alpha * sqrt(d1*d2*r);

	%% Define constraints
	% Use nuclear-norm constraint only
	%funProj = @(x,projTol,projData) projNucnorm(x,d1,d2,radius,projTol,projData);

	% Use nuclear-norm plus infinity-norm constraints
	funProj = @(x,projTol,projData) projectKappaTau(x,d1,d2,radius,alpha,projTol,projData);

	%% Recover estimate Mhat of M
	[Mhat,info] = spgSolver(funObj, funProj, zeros(d1*d2,1), options);
	Mhat = reshape(Mhat,d1,d2);
	%% From original, uncommented for speed considerations
	%[U,S,V] = svd(Mhat);
	%Mhat_debias = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; % Project onto actual rank if known

	%% Compute relative error
	%norm(Mhat_debias-M,'fro')/norm(M,'fro')
end
