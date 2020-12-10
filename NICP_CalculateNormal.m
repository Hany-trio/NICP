  function [Normal_Mat] = NICP_CalculateNormal( Point_Cloud )  

    %% parent:NICP  ���㷨��

    [ l , w ] = size(Point_Cloud);
    
    if ( w ~= 3 )
        
        Point_Cloud =  Point_Cloud';
        % ϣ�� Pointcloud ��һ�� n * 3 �ľ���    
        [ l , w ] = size(Point_Cloud);
        
    end
    
    
    PointCloud_kdtreeobj = KDTreeSearcher( Point_Cloud , 'distance' , 'euclidean' );  
    
    
    Normal_Mat = ones(l,3);
    
	for i = 1 : l   %��ȡi����Χ�ĵ㼯�Ϲ��ɵ㼯����õ㷨��

        
        [ idx, ~ ] = knnsearch( PointCloud_kdtreeobj , Point_Cloud(i,1:3) , 'dist' , 'euclidean' , 'k' , 40 );  
        
        CorrPoint = Point_Cloud( idx , 1:3 ); 
        
        [Nor_Param_A,Nor_Param_B,Nor_Param_C] = PlaneFitting( CorrPoint );
        
        Normal_Mat( i , 1:3 ) = [ Nor_Param_A , Nor_Param_B , Nor_Param_C ];
        
        
	end




  end
  
  
  
  
  
  function [A,B,C] = PlaneFitting( Points )
    
    %������ AX+BY+CZ+1 = 0;
    % �÷��̼���������� Ai + Bj +Ck ��ABC���Ƿ��߷���

    x = Points(:,1);y = Points(:,2);z = Points(:,3);

    temp_xx = x' * x;   %  sum( x .* x );
    temp_xy = x' * y;
    temp_yy = y' * y;
    temp_zx = z' * x;
    temp_zy = z' * y;
    temp_zz = z' * z;

    a = [temp_xx,temp_xy,temp_zx;
        temp_xy,temp_yy,temp_zy;
        temp_zx,temp_zy,temp_zz];
    
    b = [ -sum(x);-sum(y);-sum(z)];

    c = ((a)^-1) * b;
    
    A = c(1) ; B = c(2); C =c(3);


  end