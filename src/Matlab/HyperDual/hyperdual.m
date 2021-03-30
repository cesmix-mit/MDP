classdef hyperdual
    
    %properties (SetAccess=private, GetAccess=public)
    properties(Access=public)
        x           %#ok<*PROP>
        dx1
        dx2
        dx1x2
    end
    
    methods (Access=private) % hyperdual overloaded function templates
        
        % template for a scalar function of one (potentially vectorial) argument, i.e.      f : R^n --> R^1
        function [lhs] = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd)
            % read right-hand side
            if( ~isa(rhs,'hyperdual') )
                rhs             = hyperdual(rhs);
            end
            x               = rhs.x; 
            dx1             = rhs.dx1;
            dx2             = rhs.dx2;
            dx1x2           = rhs.dx1x2;
            [n,~,M]         = size(x);    % n is input dimension; M is SIMD page dimension
            % initialise 0th, 1st and 2nd derivative arrays
            f               = zeros(1,1,M);
            Df              = zeros(1,n,M);
            Hf              = zeros(n,n,M);
            % make the SIMD call
            if(nargin<5)
                opt_simd        = 'sequential';
            end
            switch opt_simd
                case 'sequential'
                    for k=1:M
                        f( 1,1,k)       = fun_f(  x(:,1,k) );
                        Df(1,:,k)       = fun_Df( x(:,1,k) );
                        Hf(:,:,k)       = fun_Hf( x(:,1,k) );
                    end
                case 'elemental'
                    f( 1,1,:)       = fun_f(  x(:,1,:) );
                    Df(1,:,:)       = fun_Df( x(:,1,:) );
                    Hf(:,:,:)       = fun_Hf( x(:,1,:) );
            end % switch
            % create output arrays
            m           = 1;    % the output dimension of a scalar function is "1".
            switch 'tune_1'
                case 'naive'
                    z               = zeros(m,1,M);
                    dz1             = zeros(m,1,M);
                    dz2             = zeros(m,1,M);
                    dz1z2           = zeros(m,1,M);
                    %   Unfortunately, the product below requires the
                    % for-loop since Matlab does not provide a SIMD page
                    % dimension.
                    %   This is horribly slow. Thus, there are tuned
                    % variants of this code fragment within the switch.
                    for k=1:M
                        z(    1,1,k)    = f(:,:,k);
                        dz1(  1,1,k)	= Df(1,:,k) * dx1(  :,1,k);
                        dz2(  1,1,k)  	= Df(1,:,k) * dx2(  :,1,k);
                        dz1z2(1,1,k)  	= dx1(:,1,k).' * Hf(:,:,k) * dx2(:,1,k) + Df(1,:,k) * dx1x2(:,1,k);
                    end
                case 'tune_1'   
                    %   attempt to improve the performance of the above
                    % computation when n<<M.
                    %   This is achieved by an elemental (SIMD) call over a
                    % longer array dimension, namely of M. As a result,
                    % there are less for-iterations, each of higher
                    % computational intensity.
                    z               = f;
                    dz1             = zeros(m,1,M);
                    dz2             = zeros(m,1,M);
                    dz1z2           = zeros(m,1,M);
                    for k1=1:n
                        dz1             = dz1   + Df(1,k1,:) .*   dx1(k1,1,:);
                        dz2             = dz2   + Df(1,k1,:) .*   dx2(k1,1,:);
                        dz1z2         	= dz1z2 + Df(1,k1,:) .* dx1x2(k1,1,:);
                        for k2=1:n
                            dz1z2(1,1,:)    = dz1z2(1,1,:) + dx1(k1,1,:) * Hf(k1,k2,:) * dx2(k2,1,:);
                        end
                    end
            end % switch
            % write lhs
            lhs         = hyperdual(z,dz1,dz2,dz1z2);
        end
        
        function [lhs] = fun_scalar_bivariate(fun_g,fun_DgDx,fun_DgDy,...
                            fun_DDgDxx,fun_DDgDxy,fun_DDgDyy,rhs1,rhs2,opt_simd)
            % g is still a scalar function, namely:
            %       g : R x R --> R, (x,y) |--> z=g(x,y)
            %
            % Thus we can use "fun_scalar". In order to do so, we define a
            % scalar function f.
            %       f : R^2 --> R, X |--> z=f(X):=g( X(1),X(2) )
            fun_f       = @(X)       [fun_g(      X(1,:,:) , X(2,:,:) ) ];       %#ok<NBRAK>
            fun_Df      = @(X)       [fun_DgDx(   X(1,:,:) , X(2,:,:) ) , fun_DgDy(   X(1,:,:) , X(2,:,:) ) ];
            fun_Hf    	= @(X)       [fun_DDgDxx( X(1,:,:) , X(2,:,:) ) , fun_DDgDxy( X(1,:,:) , X(2,:,:) ) ;
                                      fun_DDgDxy( X(1,:,:) , X(2,:,:) ) , fun_DDgDyy( X(1,:,:) , X(2,:,:) ) ];
            % make vectorised hyperdual
            if( ~isa(rhs1,'hyperdual') )
                rhs1        = hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2        = hyperdual(rhs2);
            end
            rhs         = hyperdual( [rhs1.x        ;rhs2.x         ],...
                                     [rhs1.dx1      ;rhs2.dx1       ],...
                                     [rhs1.dx2      ;rhs2.dx2       ],...
                                     [rhs1.dx1x2    ;rhs2.dx1x2     ] );
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
    end
    
    methods (Access=private) % private power functions
        
        function [lhs] = exponial(c,rhs)
            % c^rhs
            fun_f       = @(x)  c.^x;
            fun_Df      = @(x)  c.^x .* log(c);
            fun_Hf      = @(x)  c.^x .* log(c).^2;
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = monomial(rhs,a)
            % rhs^a
            fun_f       = @(x)  x.^a;
            fun_Df      = @(x)  a.*x.^(a-1);
            fun_Hf      = @(x)  a.*(a-1).*x.^(a-2);
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = powerial(rhs1,rhs2)
            % f(x,y):= x^y
            fun_f       = @(x,y)  x.^y;
            fun_DfDx  	= @(x,y)  y.*x.^(y-1);
            fun_DfDy  	= @(x,y)  x.^y .* log(x);
            fun_DDfDxx	= @(x,y)  y.*(y-1).*x.^(y-2);
            fun_DDfDxy	= @(x,y)  x.^(y-1).*( 1+y.*log(x) );
            fun_DDfDyy	= @(x,y)  x.^y .* log(x).^2;
            opt_simd    = 'elemental';
            lhs         = fun_scalar_bivariate( fun_f,fun_DfDx,fun_DfDy,...
                fun_DDfDxx,fun_DDfDxy,fun_DDfDyy,rhs1,rhs2,opt_simd );
        end
        
    end
    
    methods (Access=public) % class def
        
        % constructor
        function [lhs] = hyperdual(x,dx1,dx2,dx1x2)
            if(isempty(x) || size(x,1)==0)
                error('input for x must be at least one numerical value.');
            end
            size_x          = size(x);
            if(nargin<2)
                dx1             = zeros(size_x);
            else
                size_dx1        = size(dx1);
                if(norm(size_dx1-size_x,'inf')~=0)
                    error('dx1 has another size than x.');
                end
            end
            if(nargin<2)
                dx2             = zeros(size_x);
            else
                size_dx2        = size(dx2);
                if(norm(size_dx2-size_x,'inf')~=0)
                    error('dx2 has another size than x.');
                end
            end
            if(nargin<3)
                dx1x2           = zeros(size_x);
            else
                size_dx1x2     	= size(dx1x2);
                if(norm(size_dx1x2-size_x,'inf')~=0)
                    error('dx1x2 has another size than x.');
                end
            end
            % at this point it is guaranteed that the input is valid.
            lhs.x           = x;
            lhs.dx1         = dx1;
            lhs.dx2         = dx2;
            lhs.dx1x2       = dx1x2;
        end
        
        % number access
        function [varargout] = size(rhs,varargin)
            x   = rhs.x;
            [varargout{1:nargout}] = size(x, varargin{:});
        end
        
        function [res] = numel(rhs)
            res     = prod( size(rhs) ); %#ok<PSIZE>
        end
        
        function varargout = subsref(obj,s)
            switch s(1).type
                case '.'
                    if length(s) == 1
                        % Implement obj.PropertyName
                        varargout{:}    = obj.(s.subs);
                    elseif length(s) == 2 && strcmp(s(2).type,'()')
                        % Implement obj.PropertyName(indices)
                        tmp             = obj.(s(1).subs);
                        varargout{:}    = tmp(s(2).subs{:});
                    else
                        varargout       = {builtin('subsref',obj,s)};
                    end
                case '()'
                    if length(s) == 1
                        % Implement obj(indices)
                        varargout{:}    = hyperdual( obj.x(s(1).subs{:}),...
                            obj.dx1(s(1).subs{:}),obj.dx2(s(1).subs{:}),obj.dx1x2(s(1).subs{:}) );
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj(ind).PropertyName
                        % Is identical to obj.PropertyName(indices)
                        tmp             = obj.(s(2).subs);
                        varargout{1}    = tmp(s(1).subs{:});
                    elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
                        % Implement obj(indices).PropertyName(indices)
                        % Shall not exist
                        error('Multiple subindices have been used in hyperdual.');
                    else
                        % Use built-in for any other expression
                        varargout       = {builtin('subsref',obj,s)};
                    end
                case '{}'
                    error('Cell array access to hyperdual is not defined.');
                    if length(s) == 1 %#ok<UNRCH>
                        % Implement obj{indices}
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj{indices}.PropertyName
                    else
                        % Use built-in for any other expression
                        varargout = {builtin('subsref',obj,s)};
                    end
                otherwise
                    error('Not a valid indexing expression')
            end % switch
        end
        
        function [lhs] = vertcat(rhs1,rhs2)
            x       = rhs1.x;
            dx1     = rhs1.dx1;
            dx2     = rhs1.dx2;
            dx1x2   = rhs1.dx1x2;
            y       = rhs2.x;
            dy1     = rhs2.dx1;
            dy2     = rhs2.dx2;
            dy1y2   = rhs2.dx1x2;
            lhs     = hyperdual([x;y],[dx1;dy1],[dx2;dy2],[dx1x2;dy1y2]);
        end
        
        function [lhs] = horzcat(rhs1,rhs2)
            x       = rhs1.x;
            dx1     = rhs1.dx1;
            dx2     = rhs1.dx2;
            dx1x2   = rhs1.dx1x2;
            y       = rhs2.x;
            dy1     = rhs2.dx1;
            dy2     = rhs2.dx2;
            dy1y2   = rhs2.dx1x2;
            lhs     = hyperdual([x,y],[dx1,dy1],[dx2,dy2],[dx1x2,dy1y2]);
        end
        
    end
    
    % unary operators
    methods (Access=public)
        
        % sign
        function [lhs] = uplus(rhs)
          	lhs         = hyperdual(rhs.x,rhs.dx1,rhs.dx2,rhs.dx1x2);
        end
        
        function [lhs] = uminus(rhs)
            lhs         = hyperdual(-rhs.x,-rhs.dx1,-rhs.dx2,-rhs.dx1x2);
        end
        
        function [res] = sign(rhs)
            res         = sign( rhs.x );
        end
        
        function [lhs] = abs(rhs)
            if( sign(rhs)>=0 )
                lhs     = +rhs;
            else
                lhs     = -rhs;
            end
        end
        
        % transpose
        function [lhs] = ctranspose(rhs) % a'
            x           = rhs.x;
            dx1         = rhs.dx1;
            dx2         = rhs.dx2;
            dx1x2       = rhs.dx1x2;
            x           = conj( permute(x    ,[2,1,3]) );
            dx1         = conj( permute(dx1  ,[2,1,3]) );
            dx2         = conj( permute(dx2  ,[2,1,3]) );
            dx1x2       = conj( permute(dx1x2,[2,1,3]) );
            lhs         = hyperdual(x,dx1,dx2,dx1x2);
        end
        
        function [lhs] = transpose(rhs) % a.'
            x           = rhs.x;
            dx1         = rhs.dx1;
            dx2         = rhs.dx2;
            dx1x2       = rhs.dx1x2;
            x           = permute(x    ,[2,1,3]);
            dx1         = permute(dx1  ,[2,1,3]);
            dx2         = permute(dx2  ,[2,1,3]);
            dx1x2       = permute(dx1x2,[2,1,3]);
            lhs         = hyperdual(x,dx1,dx2,dx1x2);
        end
        
    end
    
    % binary operators
    methods (Access=public) 
        
        function [lhs] = plus(rhs1,rhs2)
            % f(x,y):= x+y
            %   I know it's too simple but I love brute force.
            fun_f       = @(x,y)  x+y;
            fun_DfDx  	= @(x,y)  1;
            fun_DfDy  	= @(x,y)  1;
            fun_DDfDxx	= @(x,y)  0;
            fun_DDfDxy	= @(x,y)  0;
            fun_DDfDyy	= @(x,y)  0;
            opt_simd    = 'elemental';
            lhs         = fun_scalar_bivariate( fun_f,fun_DfDx,fun_DfDy,...
                fun_DDfDxx,fun_DDfDxy,fun_DDfDyy,rhs1,rhs2,opt_simd );
        end
        
        function [lhs] = minus(rhs1,rhs2)
            lhs     = rhs1 + (-rhs2);   % exploit unary minus
        end
        
        function [lhs] = mtimes(rhs1,rhs2) % * 
            % f(x,y):= x*y
            if( ~isa(rhs1,'hyperdual') )
                if( isa(rhs1,'double') && size(rhs1,1)==1 && size(rhs1,2)==1 )
                    rhs1            = repmat(rhs1,size(rhs2));
                else
                    rhs1            = repmat( rhs1(:,:,1), [1,1,size(rhs2,3)] );
                end
             	rhs1           	= hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2          	= hyperdual(rhs2);
            end
            % read
            X       = rhs1.x;
            dX1    	= rhs1.dx1;
            dX2     = rhs1.dx2;
            dX1X2   = rhs1.dx1x2;
            Y       = rhs2.x;
            dY1    	= rhs2.dx1;
            dY2     = rhs2.dx2;
            dY1Y2   = rhs2.dx1x2;
            % write:
            %   We use the direct formula and take care that the expressions
            %   work well also when X and Y are matrices.
            M       = size(X,3);
            Z       = zeros(size(X,1),size(Y,2),M);
            dZ1     = zeros(size(X,1),size(Y,2),M);
            dZ2     = zeros(size(X,1),size(Y,2),M);
            dZ1Z2   = zeros(size(X,1),size(Y,2),M);
            for k=1:M   % simd direction
                Z(     :,:,k) 	= X(:,:,k) * Y(  :,:,k);
                dZ1(   :,:,k) 	= X(:,:,k) * dY1(:,:,k) + dX1(:,:,k) * Y(:,:,k);
                dZ2(   :,:,k) 	= X(:,:,k) * dY2(:,:,k) + dX2(:,:,k) * Y(:,:,k);
                dZ1Z2( :,:,k)   = X(:,:,k) * dY1Y2(:,:,k) + dX1(:,:,k) * dY2(:,:,k) + dX2(:,:,k) * dY1(:,:,k) + dX1X2(:,:,k) * Y(:,:,k);
            end
            lhs     = hyperdual(Z,dZ1,dZ2,dZ1Z2);
        end
        
        function [lhs] = times(rhs1,rhs2) % .*
            % f(x,y):= x.*y
            if( ~isa(rhs1,'hyperdual') )
                rhs1           	= hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2          	= hyperdual(rhs2);
            end
            % read
            X       = rhs1.x;
            dX1    	= rhs1.dx1;
            dX2     = rhs1.dx2;
            dX1X2   = rhs1.dx1x2;
            Y       = rhs2.x;
            dY1    	= rhs2.dx1;
            dY2     = rhs2.dx2;
            dY1Y2   = rhs2.dx1x2;
            lhs     = hyperdual(...
                                X.*Y,...
                                X.*dY1+dX1.*Y,...
                                X.*dY2+dX2.*Y,...
                                dX1.*dY2+dX2.*dY1+X.*dY1Y2+dX1X2.*Y...
                                    );
        end
        
        function [lhs] = rdivide(rhs1,rhs2)% ./
            % f(x,y):= x./y
            %   We could use power and times but this would both
            %       - be more computationally expensive
            %       - introduce more round-off
            %
            fun_f       = @(x,y)  x./y;  % x and y are elemental
            fun_DfDx  	= @(x,y)  1./y;
            fun_DfDy  	= @(x,y) -x./y.^2;
            fun_DDfDxx	= @(x,y)  0*x;
            fun_DDfDxy	= @(x,y) -1./y.^2;
            fun_DDfDyy	= @(x,y)  2*x./y.^3;
            opt_simd    = 'elemental';
            lhs         = fun_scalar_bivariate( fun_f,fun_DfDx,fun_DfDy,...
                fun_DDfDxx,fun_DDfDxy,fun_DDfDyy,rhs1,rhs2,opt_simd );
        end
        
        function [lhs] = ldivide(rhs1,rhs2)% .\
            [lhs] = rdivide(rhs2,rhs1);
        end
        
        function [lhs] = mrdivide(rhs1,rhs2) %#ok<STOUT,MANU,INUSD>
            % f(x,y):= x/y
            %   not defined on purpose!
            %   Reason: We only attempt to overload functions of scalar in-
            %   and output.
            error('The operator MRDIVIDE is not defined for hyperdual.\n');
        end
        
        function [lhs] = mldivide(rhs1,rhs2) %#ok<STOUT,MANU,INUSD>
            % f(x,y):= x\y
            %   not defined on purpose!
            %   Reason: We only attempt to overload functions of scalar in-
            %   and output.
            error('The operator MLDIVIDE is yet not defined for hyperdual.\n');
        end
        
        function [lhs] = mpower(rhs1,rhs2) %#ok<STOUT,MANU,INUSD>
            % f(x,y):= x^y
            %   not defined on purpose!
            %   Reason: We only attempt to overload functions of scalar in-
            %   and output.
            error('The operator MPOWER is yet not defined for hyperdual.\n');
        end
        
        function [lhs] = lt(rhs1,rhs2)% <
            if( ~isa(rhs1,'hyperdual') )
                rhs1           	= hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2          	= hyperdual(rhs2);
            end
            lhs     = (rhs1.x) < (rhs2.x);
        end
        
        function [lhs] = gt(rhs1,rhs2)% >
            if( ~isa(rhs1,'hyperdual') )
                rhs1           	= hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2          	= hyperdual(rhs2);
            end
            lhs     = (rhs1.x) > (rhs2.x);
        end
        
        function [lhs] = le(rhs1,rhs2)% <=
            if( ~isa(rhs1,'hyperdual') )
                rhs1           	= hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2          	= hyperdual(rhs2);
            end
            lhs     = (rhs1.x) <= (rhs2.x);
        end
        
        function [lhs] = ge(rhs1,rhs2)% >=
            if( ~isa(rhs1,'hyperdual') )
                rhs1           	= hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2          	= hyperdual(rhs2);
            end
            lhs     = (rhs1.x) >= (rhs2.x);
        end
        
        function [lhs] = ne(rhs1,rhs2)% ~=
            if( ~isa(rhs1,'hyperdual') )
                rhs1           	= hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2          	= hyperdual(rhs2);
            end
            lhs     = (rhs1.x) ~= (rhs2.x);
        end
        
        function [lhs] = eq(rhs1,rhs2)% ==
            if( ~isa(rhs1,'hyperdual') )
                rhs1           	= hyperdual(rhs1);
            end
            if( ~isa(rhs2,'hyperdual') )
                rhs2          	= hyperdual(rhs2);
            end
            lhs     = (rhs1.x) == (rhs2.x);
        end
                
    end
    
    % natural functions
    methods (Access=public) % trigonometric functions
        
        % sin,cos,tan
        function [lhs] = sin(rhs)
            fun_f       = @(x)  sin(x);
            fun_Df      = @(x)  cos(x);
            fun_Hf      = @(x) -sin(x);
            opt_simd    = 'elemental';  % the above functions are well-defined for elemental call of nD arrays.
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = cos(rhs)
            fun_f       = @(x)  cos(x);
            fun_Df      = @(x) -sin(x);
            fun_Hf      = @(x) -cos(x);
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = tan(rhs)
            fun_f       = @(x)  tan(x);
            fun_Df      = @(x)  sec(x).^2;
            fun_Hf      = @(x)  2*tan(x).*sec(x).^2;
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        % arcsin,arccos,arctan (area, i.e. inverse; it occurs Matlab calls them "asin,acos,atan")
        function [lhs] = asin(rhs)
            fun_f       = @(x)  asin(x);
            fun_Df      = @(x)  1./(1-x.^2).^(1/2);
            fun_Hf      = @(x)  x./(1-x.^2).^(3/2);
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = acos(rhs)
            fun_f       = @(x)  acos(x);
            fun_Df      = @(x) -1./(1-x.^2).^(1/2);
            fun_Hf      = @(x) -x./(1-x.^2).^(3/2);
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = atan(rhs)
            fun_f       = @(x)  atan(x);
            fun_Df      = @(x)  1./(1+x.^2);
            fun_Hf      = @(x) -2*x./(1+x.^2).^2;
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        % sinh,cosh,tanh (hyperbolicus)
        function [lhs] = sinh(rhs)
            fun_f       = @(x)  sinh(x);
            fun_Df      = @(x)  cosh(x);
            fun_Hf      = @(x)  sinh(x);
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = cosh(rhs)
            fun_f       = @(x)  cosh(x);
            fun_Df      = @(x)  sinh(x);
            fun_Hf      = @(x)  cosh(x);
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = tanh(rhs)
            fun_f       = @(x)  tanh(x);
            fun_Df      = @(x)  sech(x).^2;
            fun_Hf      = @(x) -2*tanh(x).*sech(x).^2;
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        % asinh,acosh,atanh (area hyperbolicus, i.e. inverse hyperbolicus)
        function [lhs] = asinh(rhs)
            fun_f       = @(x)  asinh(x);
            fun_Df      = @(x)  1./(x.^2 + 1);
            fun_Hf      = @(x) -x./(x.^2 + 1).^(3/2);
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = acosh(rhs)
            fun_f       = @(x)  acosh(x);
            fun_Df      = @(x)  1./( (x-1).^(1/2) .* (x+1).^(1/2) );
            fun_Hf      = @(x) -x./( (x-1).^(3/2) .* (x+1).^(3/2) );
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = atanh(rhs)
            fun_f       = @(x)  acosh(x);
            fun_Df      = @(x)  1./( 1-x.^2 );
            fun_Hf      = @(x)  2*x./( 1-x.^2 ).^2;
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
    end
    
    methods (Access=public) % power functions
        
        % exp, log (i.e. lagartihmus naturalis)
        function [lhs] = exp(rhs)
            fun_f       = @(x)  exp(x);
            fun_Df      = @(x)  exp(x);
            fun_Hf      = @(x)  exp(x);
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        function [lhs] = log(rhs)
            fun_f       = @(x)  log(x);
            fun_Df      = @(x)  1./x;
            fun_Hf      = @(x) -1./x.^2;    % ( ashamed to have used W-alpha for that :) )
            opt_simd    = 'elemental';
            lhs         = fun_scalar(fun_f,fun_Df,fun_Hf,rhs,opt_simd);
        end
        
        % power
        % For variables x,y and constants c,a we have to separate the
        % operations 
        %       1) c^y      "exponial"
        %       2) x^a      "monomial"
        %       3) x^y      "powerial" :)
        % because of the following reasons:
        %   A) If basis==0 or power==0 then d/dx is indeterminate for (1&3) whereas well-defined for (2).
        %   B) If basis<0 and power not integer then basis^power is indeterminate.
        %   C) Since (B) the expressions (1&3) do only make sense for x>0
        function [lhs] = power(rhs1,rhs2)
            b_hd1           = isa(rhs1,'hyperdual');
            b_hd2           = isa(rhs2,'hyperdual');
            b_exponial      = b_hd2 && (~b_hd1);
            b_monomial      = b_hd1 && (~b_hd2);
            b_powerial      = b_hd1 && b_hd2;
            if( b_exponial )
                lhs             = exponial(rhs1,rhs2);
            end
            if( b_monomial )
                lhs             = monomial(rhs1,rhs2);
            end
            if( b_powerial )
                lhs             = powerial(rhs1,rhs2);
            end
        end
        
    end
end