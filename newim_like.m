% newim_like : chooses the appropriate newim function depending on the input
% Decides whether the standard dipimage function or the cuda version is used, depending on the state of the global use_newim_cuda variable, 
% which can be set via the function set_newim_cuda()

%***************************************************************************
%   Copyright (C) 2008-2009 by Rainer Heintzmann                          *
%   heintzmann@gmail.com                                                  *
%                                                                         *
%   This program is free software; you can redistribute it and/or modify  *
%   it under the terms of the GNU General Public License as published by  *
%   the Free Software Foundation; Version 2 of the License.               *
%                                                                         *
%   This program is distributed in the hope that it will be useful,       *
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
%   GNU General Public License for more details.                          *
%                                                                         *
%   You should have received a copy of the GNU General Public License     *
%   along with this program; if not, write to the                         *
%   Free Software Foundation, Inc.,                                       *
%   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
%**************************************************************************
%

function res=newim_like(varargin)
global diphandle_newim;
if isa(varargin{1},'cuda')
    tmp=cuda();  % this serves to call the cuda version
    res= newim_cuda2(tmp,varargin{2:end});        
else
    res=feval(diphandle_newim,varargin{2:end});  % call the standart matlab zeros function
end
