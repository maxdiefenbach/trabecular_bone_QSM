classdef MyHandle < matlab.mixin.Copyable & dynamicprops
% base class to handle property io and version trace back
    
    properties
        Version
        fileID
    end
    
    methods
        
        function this = MyHandle
            if isempty(this.Version)
                this.set_Version;
            end
        end
        
        
        function Version = set_Version(this)
        % set Version Struct with git SHA1, matlab version, and datetime char
        % user and computer name (host)
            fpath = fileparts(mfilename('fullpath'));
            Version.BMRR = get_GitInfo(fpath);
            Version.matlab = version;
            Version.datetime = char(datetime);
            try
                [~, USER] = system('echo $USER');
                Version.user = USER;
                [~, HOST] = system('echo $HOST');
                Version.host = HOST;
            catch
                fprintf('Implement me for Windows.')
            end
            this.Version = Version;
        end
        
        
        function set_fileID(this, fileID)
            if nargin == 2
                this.fileID = fileID;
            end
        end

        
        function load_propertyFromFile(this, propertyName, filename)
        % load_propertyFromFile(this, propertyName, filename)
        % 
        % propertyName, char
        % filename, char, possible path to file or dir
        % if isdir(filename) look for file with the name <filename>/<fileID>_<propertyName>.mat
            
            if nargin < 3
                fpath = fileparts(this.filename);
                filename = fullfile(fpath, [this.fileID '_' propertyName '.mat']);
            elseif isdir(filename)
                [fpath, fname] = fileparts(filename);
                filename = fullfile(fpath, fname, [this.fileID '_' propertyName '.mat'])
            end

            tmp = load(filename);
            this.(propertyName) = tmp.(propertyName)
            
            S = dbstack;
            fprintf('loaded %s from %s (%s:%d)\n', propertyName, GetFullPath(filename), mfilename, S(1).line+1);
        end

        
        function filename = save_property(this, propertyName, filename, version, overwrite)
        % save_property(this, propertyName, filename, version) 
        % 
        % save a property with propertyName to matfile filename with mat-file version version 
        % adds Version struct to property struct
        % 
        % propertyName char, possible with dot notation to indicate several levels,
        %              e.g. 'Struct.field1.field2' 
        % filename char, possible path to file or directory if director does not exist, it
        %          gets created 
        % version char, default '-v7.3'
            
            assert(~isempty(this.(propertyName)));
            if nargin < 5
                overwrite = 1;
            end
            if nargin < 4
                version = '-v7.3';      % default mat file version
            end
            propertyNameCell = strsplit(propertyName, '.');
            if nargin < 3
                filename = [this.fileID '_' propertyNameCell{1} '.mat'];
            end
            [~, fname, fext] = fileparts(filename);
            if isdir(filename)          % if dir exists 
                filename = fullfile(filename, [this.fileID '_' propertyNameCell{1} '.mat']);
            elseif ~isdir(filename) & isempty(fext) % if dir does not exist
                mkdir(filename);
                filename = fullfile(filename, [this.fileID '_' propertyNameCell{1} '.mat']);
            end
            if exist(filename, 'file')
                S = dbstack;
                fprintf('File %s already exists. (%s:%d)\n', fname, mfilename, S(1).line+1)
                if overwrite
                    S = dbstack;
                    fprintf('Overwrite. (%s:%d)\n', mfilename, S(1).line+1)
                else
                    error('Protect.\n')
                end
            end
            
            nLevels = numel(propertyNameCell);
            cmd = 'Struct';
            propertyChain = '';
            for iLevel = 1:nLevels
                propertyChain = [propertyChain sprintf('.(propertyNameCell{%d})', iLevel)];
            end
            cmd = [cmd propertyChain ' = this' propertyChain ';'];
            eval(cmd);
            
            assert(~isempty(Struct))
            Struct.Version = this.Version;
            save(filename, '-struct', 'Struct', version); % puts Version next to Struct in the same file
            
            S = dbstack;
            fprintf('wrote %s (%s:%d)\n', GetFullPath(filename), mfilename, S(1).line+1);
        end
        
        
        function prop = get_property(this, propName)
            propNameCell = strsplit(propName, '.');
            nLevels = numel(propNameCell);
            propChain = '';
            for iLevel = 1:nLevels
                propChain = [propChain sprintf('.(propNameCell{%d})', iLevel)];
            end
            cmd = ['prop = this' propChain ';'];
            eval(cmd);
            assert(~isempty(prop));
        end

    end % methods

end % classdef