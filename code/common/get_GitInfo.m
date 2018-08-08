function GitInfo = get_GitInfo(path)
    
    if exist(path, 'file') & ~isdir(path) % if path is a filename
        path = fileparts(path);
    end

    if ~exist(fullfile(path, '.git/'), 'dir') % move to parent directory
        uppath = fileparts(path);
        if strcmp(uppath, '/')          % prevent infinite recursion
            GitInfo = struct();
            return
        end
        GitInfo = get_GitInfo(uppath);
        return
    end
    
    head = fileread(fullfile(path, '.git/HEAD'));
    tmp = strsplit(head, '/');
    branch = char(strrep(tmp(end), char(10), '')); % remove newline character ^Z

    SHA1file = fileread(fullfile(path, '.git/refs/heads', branch));
    commit = char(strrep(SHA1file, char(10), ''));
    
    remote = '';
    url = '';
    config = fullfile(path, '.git/config');
    fid = fopen(config);
    nLines = length(config);
    for iL = 1:nLines
        line = fgetl(fid);
        if isequal(line, sprintf('[branch "%s"]', branch))
            remoteline = fgetl(fid);
            tmp = strsplit(remoteline);
            remote = char(tmp(end));
        end
    end
    fclose(fid);
    
    fid = fopen(config);
    for iL = 1:nLines
        line = fgetl(fid);
        if isequal(line, sprintf('[remote "%s"]', remote))
            urlline = fgetl(fid);
            tmp = strsplit(urlline);
            url = char(tmp(end));
        end
    end
    fclose(fid);
    
    GitInfo.branch = branch;
    GitInfo.commit = commit;
    GitInfo.remote = remote;
    GitInfo.url = url;

end