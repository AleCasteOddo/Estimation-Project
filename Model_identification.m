function [identification varargout] = Model_identification(simulation_data,sample_time,model_fun,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation and Learning in Aerospace Project A.Y. 24-25 
% Function to identify by grey_box ID

% Authors:  Alessandro Castelanelli Oddo (alessandro.castelanelli@polimi.it)
%            (@polimi.it)                     
%            (@polimi.it)                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data in-take
data.ax = simulation_data.ax;
data.q = simulation_data.q;
input = simulation_data.Mtot;

% Measured acceleration and pitch rate
output = [data.q data.ax];
% Data ordering and bring those in frequency domain
sim_data = iddata(output, input, sample_time);
data_fd = fft(sim_data); % output of the simulation in the frequency domain

% Initial guess for the identification
guess = zeros(6,1);
sys_init = idgrey(model_fun, guess, 'c');

% Actual Model Identification
identification = struct;
estimated_model = greyest(data_fd, sys_init);
identification.parameters = estimated_model.Report.Parameters.ParVector;
identification.fit = estimated_model.Report.Fit.FitPercent;
identification.covariance = getcov(estimated_model);
identification.matrix={estimated_model.A; estimated_model.B; estimated_model.C; estimated_model.D};
identification.estimated_model=estimated_model;

% Data output
if nargout>=2
    real_parameters = varargin{1};
    varargout{1} = (identification.parameters-real_parameters) ./ real_parameters * 100; %Estimation Error
end
if nargout>=3
    varargout{2} = [input output];
end
if nargout==4
    A=identification.matrix{1};
    B=identification.matrix{2};
    C = [1 0 0 ; identification.matrix{3}(1,:) ; 0 0 1 ; identification.matrix{3}(2,:)];
    D = [0; 0 ; identification.matrix{4}];
    sim_est.ax = simulation_data.ax.Data;
    sim_est.q = simulation_data.q.Data;
    varargout{3} = sim_est;
end

end
