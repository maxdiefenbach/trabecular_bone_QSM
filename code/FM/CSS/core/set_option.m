function val = set_option(Options, fieldName, defaultVal)
    
    if nargin < 3
        defaultVal = [];
    end
    
    if isfield(Options, fieldName) & ~isempty(Options.(fieldName))
        val = Options.(fieldName);

%   In codegen 
% ??? The function 'isprop' is not supported for standalone code generation. See the documentation for
%         coder.extrinsic to learn how you can use this function in simulation.
% elseif isprop(Options, fieldName) & ~isempty(Options.(fieldName))
%     val = Options.(fieldName);
    
    else
        val = defaultVal;
%     if isscalar(defaultVal)
%         fprintf('set %s to %s\n', fieldName, num2str(defaultVal))
%     else
%         if isempty(defaultVal)
%             fprintf('set %s to empty\n', fieldName);
%         else
%             fprintf('set %s to shape [%s], el(1, 1) = %s\n', fieldName, num2str(size(defaultVal)), num2str(defaultVal(1, 1)));
%         end
% end

    end
