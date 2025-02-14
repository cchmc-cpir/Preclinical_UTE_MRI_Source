classdef ReconModel
	properties
		verbose;
		unique_string;
        system;
        deapodize = true;
        %PJN-Add for Fermi filter - default to no
        fermifilter = false;
        crop = true;
	end
	
	methods
		function obj = ReconModel(system_model, verbosity)
			obj.verbose = verbosity;
			obj.system = system_model;
        end
    end
    
    methods (Abstract)
		% Reconstructs an image volume using the given data
		reconVol = reconstruct(obj,data, traj);
	end
end
