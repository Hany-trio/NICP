%% NICP algor
    %20201209 Hany Njtech
 
function [T] = HCode_NICP(Refer_cloud,Match_cloud,Refer_Normal,Match_Normal)
    
    % 输入 ： 点云 1 2 文件； 返回： T 变换矩阵
    % Refer_cloud 参考点云 N*3 不做变换操作
    % Match_cloud 匹配点云 N*3 做变换操作
    % Refer_Normal 参考点云法线 N*3
    % Match_Normal 匹配点云法线 N*3
    
    
    %% 计算初始平移系数
    
    [ Refer_cloud_l , ~ ] = size(Refer_cloud);  %读取点云大小
    [ Match_cloud_l , ~ ] = size(Match_cloud);
    
    CoG_Refer_cloud = sum(Refer_cloud) / Refer_cloud_l;  %计算点云重心
    CoG_Match_cloud = sum(Match_cloud) / Match_cloud_l;    
    
    Initial_Translation_Param = CoG_Refer_cloud - CoG_Match_cloud;  %初始平移量
    Initial_Rotation_Param = diag(ones(1,3));                       %初始旋转量
    
   
  %% 先进行初始变换    
    
    R = Initial_Rotation_Param;   % 对平移以及旋转进行初始化
    t = [0,0,0]'; 
    T(1:3,1:3) = R;
    T(1:3,4) = t;
    T(4,1:4) = [0,0,0,1];
    
    
   	%% 计算邻近关系 取欧式距离最小值作为对应点 
    
	Refer_kdtreeobj = KDTreeSearcher(Refer_cloud, 'distance' , 'euclidean' );    %建立 Refercloud的kdtree    
    
    
    %% 迭代计算
    
    correspondence = zeros(Match_cloud_l,1);    
    
    P_number = 0.3;
    
    P =[diag([1,1,1]),zeros(3);
        zeros(3),diag([P_number,P_number,P_number])];
    
    
    while(true)   %迭代计算开始
        
    H = zeros(6);  % 定义 Hessian 矩阵
    b = zeros(6,1); % 定义 J * e 
    
        for i = 1 : Match_cloud_l   %确定点对应关系

            [ idx, ~ ] = knnsearch( Refer_kdtreeobj , Match_cloud(i,1:3) , 'dist' , 'euclidean' , 'k' , 1); 
            correspondence(i,1) = idx; 

        end
        
        for i = 1 : Match_cloud_l
            
            J = [  eye(3)   ,  -ToAntisymmetric_Mat( Match_cloud(i,1:3)' );
                   zeros(3) ,  -ToAntisymmetric_Mat( Refer_Normal( correspondence(i,1) , 1:3 )' ) * ToAntisymmetric_Mat( Match_Normal(i,1:3)') ]; 
            
            H = H + ( J' * P * J ); % Hessian
            
            b = b + J' * P * [ ( Match_cloud(i,1:3)' - Refer_cloud( correspondence(i,1) , 1:3)' ); ( ToAntisymmetric_Mat( Refer_Normal( correspondence(i,1) ,1:3)' ) * Match_Normal(i,1:3)')];
            

        end      
        
        delt = H \ -b;         
        
        t_ = delt(1:3,1);
        
        fa = delt(4:6,1);
        
        norm_delt = norm(fa);
        a = fa / norm_delt;
        
        J_ = ( sin( norm_delt ) / norm_delt ) * eye(3) + ( 1 - sin(norm_delt) / norm_delt )*( a * a' )+(( 1 - cos(norm_delt) ) / norm_delt ) * ToAntisymmetric_Mat(a);
        deltR = cos( norm_delt ) * eye(3) + ( 1 - cos( norm_delt ) ) * ( a * a' ) + sin( norm_delt ) * ToAntisymmetric_Mat(a);
        
        
        deltT = [deltR     ,  J_ * t_;
                 zeros(1,3) , 1      ];

        T = deltT ;

        
        Match_cloud = Trans( Match_cloud , T );
        Match_Normal = NICP_CalculateNormal( Match_cloud );
        Match_Normal = NormalGeneral(Match_Normal);
        
        close all
        pcshow(Match_cloud , [0,1,0] , 'MarkerSize' , 20);
        hold on 
        pcshow(Refer_cloud , [1,0,0] , 'MarkerSize' , 20);
        
    end
    
    


end

function [cloud] = Trans(cloud,T)

    [l , w] = size(cloud);
    
    if ( w ~= 3 )
        
        cloud =cloud';
        [l , ~] = size(cloud);
        
    end
    
    
    cloud(:,4) = ones(l,1);
    
    cloud = ( T * cloud' )';
    
    cloud = cloud( :, 1:3 );



end




function [Mat] = ToAntisymmetric_Mat(vector)
    
    [ w , ~ ] = size(vector);
    
    if(w ~= 3)
        
        vector = vector';
        
    end
       
    a1 = vector(1,1);
    a2 = vector(2,1);
    a3 = vector(3,1);
    
    
    Mat = [ 0  , -a3 ,  a2;
            a3 ,  0  , -a1;
           -a2 ,  a1 ,  0];

end


 function [outNormal] = NormalGeneral(inputNormal)
    % 让法向量标准化
    
       
    [ l , w ] = size(inputNormal);
    
    if( w ~= 3 )
        
        inputNormal = inputNormal';
        [ l , ~ ] = size(inputNormal);
        
    end
    
    outNormal = zeros(l,3);
    
    
    for i = 1 : l
        
        
        norm_Normal = norm( inputNormal( i , 1:3 ) );
        outNormal(i,1:3) = inputNormal(i,1:3) / norm_Normal;
        
    end
    
    
    
    end
