function theta = get_angular_coordinate(x,y)
    if numel(x)==0
        theta = [];
    else
        for i = 1:numel(x)
            if x(i)==0
                theta(i) = pi/2;
            else
                if x(i)>0
                    % theta(i) = atan(y(i)/x(i))-pi/2;
                    theta(i) = atan(y(i)/x(i));
                else
                    % theta(i) = atan(y(i)/x(i))+pi/2;
                    theta(i) = atan(y(i)/x(i))+pi;
                end
            end
        end
    end
end