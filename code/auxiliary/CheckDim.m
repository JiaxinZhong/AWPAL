% Check the dimensions of input arguments
function res = CheckDim(type, varargin)
    arg_num = nargin - 1;

    arg_name = cell(arg_num,1 );
    for i = 1:arg_num
        arg_name{i} = inputname(i+1);
    end
    
    switch type
        case 'preceding'
            if (numel(varargin{1}) == 1 || numel(varargin{2}) == 1)
                return
            end
            for i = 1:numel(size(varargin{1}))-1
                if (size(varargin{2}, i) > 1)
                    for j = 1:2
                        fprintf("Size of variable '%s' is %s.\n", ...
                            arg_name{j}, ...
                            strjoin(erase(cellstr(num2str(size(varargin{j}).')), ' '), 'x'));
                    end
                    error("Dimensions of '%s' should be preceding to those of '%s'!\n", ...
                        arg_name{1}, arg_name{2})
                end
            end
        otherwise
            error("not supported type!\n")
    end

%     switch arg_num
%         case 1
%         case 2
%         case 3
%         case 4
%         otherwise
%             error("Not supported number of arguments!\n");
%     end
    res = true;
end
