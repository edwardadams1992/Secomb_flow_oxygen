# Secomb_flow_oxygen
%%<br />
MATLAB implentation of flow, oxygen estimations in microvascular networks

%%<br />
MATLAB version of Secombs software for estimating flow boundary conditions,flow, and oxygen diffusion in microvascular networks.
See Secomb's original versions here:
https://github.com/secomb/NetFlowV2
https://github.com/secomb/FlowEstimateV1
https://github.com/secomb/GreensV4


%%<br />
A test network has been provided: 'test_network.mat'.
To run estimates of this network use: 'test_run.m'
To run your own networks, the structure must be defined as per 'test_network.mat'  
%%<br />
Parameters for blood, oxygen diffusion can be changed by editing constants in: 'blood.m','bloodconc.m','bloodconcp.m','convect.m','secomb_oxygen.m'
Flow parameters can be found in 'flow_estimate_bc.m', 'netflow.m', and 'secomb_oxygen.m'. 
For further details on the theory and how the alogrithm is implemented, see Secomb's work for details:
flow_estimate_bc.m - https://github.com/secomb/FlowEstimateV1
netflow.m - https://github.com/secomb/NetFlowV2
secomb_oxygen.m - https://github.com/secomb/GreensV4
%% <br />
The function 'tubeplot' is used to view flow in the networks. See: 
https://uk.mathworks.com/matlabcentral/fileexchange/5562-tubeplot

%% <br />
Matlab version used:
2022b

%%<br />
Feel free to use this code at your own risk. I welcome any feedback or acknowledgement if you find it useful.

