classdef Deque
    properties
        Data
    end

    methods
        function obj = Deque()
            obj.Data = {};
        end

        function obj = append(obj, value)
            obj.Data{end + 1} = value;
        end

        function [obj, value] = popleft(obj)
            if isempty(obj.Data)
                error('Queue is empty');
            end
            value = obj.Data{1};
            obj.Data(1) = [];
        end

        function status = is_empty(obj)
            status = isempty(obj.Data);
        end
    end

end
