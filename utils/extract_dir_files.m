function [Dir] = extract_dir_files(Dir, extractFiles,ignore_dir)
%EXTRACT_DIR extracts the directory structure of a path given by Dir.
%extractFiles = 1: only extract the files. if 0, extract the directories.
%ignore_dir (optional): flag to ignore a certain directory name. for
%example, 'all_structs'
if(nargin == 2)
    ignore_dir = [];
end
Dir = dir(Dir);
Dir = Dir(arrayfun(@(x) x.name(1), Dir) ~= '.');
if(extractFiles == 1)
    Dir = Dir(arrayfun(@(x) (x.isdir), Dir) ~=1);
else
    Dir = Dir(arrayfun(@(x) (x.isdir), Dir) ==1);
end
if(~isempty(ignore_dir))
   Dir = Dir(arrayfun(@(x) (contains(x.name,ignore_dir)), Dir) ~=1); 
end
end

