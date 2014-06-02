function [a, b, c] = getRKMethod(type, order, quadrule)
    if nargin == 2
        quadrule = 'radau';
    end
    
    switch lower(type)
        case 'dirk'
            a = getDIRKMethod(order, quadrule);        
        case 'sdirk'
            a = getSDIRKMethod(order, quadrule);
        case 'firk'
            a = getFIRKMethod(order, quadrule);
        case 'rirk'
            a = getRIRKMethod(order, quadrule);
        case 'sirk'
            a = getSIRKMethod(order, quadrule);
    end
    
    b = a(end, :);
    c = sum(a, 2);
end

function a = getDIRKMethod(order, quadrule)
    %% Diagonally-Implicit Runge-Kutta Methods
    switch quadrule
        case 'radau'
            switch order
                case 1
                    a = 1;
                case 2
                    a = zeros(2);
                    
                    a(1,1) = 1/3;
                    
                    a(2,1) = 3/4;
                    a(2,2) = 1/4;
                case 3
                    a = zeros(3);
                    
                    a(1,1) = 1/6;
                    
                    a(2,1) = 1/4;
                    a(2,2) = 1/4;
                    
                    a(3,1) = 1/6;
                    a(3,2) = -1/6;
                    a(3,3) = 1/3;
                    
                    a(4,3) = 3/4;
                    a(4,4) = 1/4;
                otherwise
                    error('DIRK-RadauIIA order must be <=3');
            end
        otherwise
            error('No DIRK methods implementated for quadrature rule %s', quadrule);
    end
end

function a = getSDIRKMethod(order, quadrule)
    %% Singuly Diagonally-Implicit Runge-Kutta Methods
    error('SDIRK methods not implemented');
end

function a = getFIRKMethod(order, quadrule)
    %% Fully-Implicit Runge-Kutta Methods
    switch quadrule
        case 'radau'
            switch order
                case 1
                    a = 1;
                case 3
                    a = zeros(2);
                    
                    a(1,1) = 5/12;
                    a(1,2) = -1/12;
                    
                    a(2,1) = 3/4;
                    a(2,2) = 1/4;
                case 5
                    a = zeros(3);
                    
                    a(1,1) = (88 - 7*sqrt(6))/360;
                    a(1,2) = (296 - 169*sqrt(6))/1800;
                    a(1,3) = (-2 + 3*sqrt(6))/225;
                    
                    a(2,1) = (296 + 169*sqrt(6))/1800;
                    a(2,2) = (88 + 7*sqrt(6))/360;
                    a(2,3) = (-2 - 3*sqrt(6))/225;
                    
                    a(3,1) = (16 - sqrt(6))/36;
                    a(3,2) = (16 + sqrt(6))/36;
                    a(3,3) = 1/9;
                otherwise
                    error('FIRK-RadauIIA order must be odd and <=5');
            end
        case 'lobatto'
            switch order
                case 2
                    a = zeros(2);
                    
                    a(1,1) = 1/2;
                    a(1,2) = -1/2;
                    
                    a(2,1) = 1/2;
                    a(2,2) = 1/2;
                case 4
                    a = zeros(3);
                    
                    a(1,1) = 1/6;
                    a(1,2) = -1/3;
                    a(1,3) = 1/6;
                    
                    a(2,1) = 1/6;
                    a(2,2) = 5/12;
                    a(2,3) = -1/12;
                    
                    a(3,1) = 1/6;
                    a(3,2) = 2/3;
                    a(3,3) = 1/6;
                case 6
                    a = zeros(4);
                    
                    a(1,1) = 1/12;
                    a(1,2) = -sqrt(5)/12;
                    a(1,3) = sqrt(5)/12;
                    a(1,4) = -1/12;
                    
                    a(2,1) = 1/12;
                    a(2,2) = 1/4;
                    a(2,3) = (10 - 7*sqrt(5))/60;
                    a(2,4) = sqrt(5)/60;
                    
                    a(3,1) = 1/12;
                    a(3,2) = (10 + 7*sqrt(5))/60;
                    a(3,3) = 1/4;
                    a(3,4) = -sqrt(5)/60;
                    
                    a(4,1) = 1/12;
                    a(4,2) = 5/12;
                    a(4,3) = 5/12;
                    a(4,4) = 1/12;
                otherwise
                    error('FIRK-LobattoIIIC order be even and <=6');
            end
        otherwise
            error('No FIRK methods implemented for quadrature rule %s', quadrule);
    end
end

function a = getRIRKMethod(order, quadrule)
    %% Real-Implicit Runge-Kutta Methods
    error('RIRK methods not implemented');
end

function a = getSIRKMethod(order, quadrule)
    %% Singly-Implicit Runge-Kutta Methods
    error('SIRK methods not implemented');
end