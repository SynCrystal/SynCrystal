

function src = syn_src(directory)
	if (nargin<1)
		directory = './scatnet/data/crystal/'; % PUT DEFAULT DIRECTORY HERE
    end
	src = create_src(directory, @uiuc_extract_objects_fun);
end

function [objects, classes] = uiuc_extract_objects_fun(file)
	objects.u1 = [1, 1];
	objects.u2 = [64, 64];
	path_str = fileparts(file);
	path_parts = regexp(path_str, filesep, 'split');
	classes = {path_parts{end}};
end
