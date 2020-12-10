%% Main 2 the second Main program
    % just solve the rotation param
    % the file just as a test

    clc
    clear
    
%% ��ȡ��������

    cloud1 = pcread("rabbit_Segment_202011142150.pcd");
%     cloud2 = pcread("rabbit_Segment_202011142150_onlyRO.pcd");   
    cloud2 = pcread("rabbit_Segment_202011142150.pcd");

    cloud1 = cloud1.Location;
    cloud2 = cloud2.Location;
   
    
    R = [0.9814  ,  0.0093  ,  0.1917;
         0.0740  ,  0.9032  , -0.4228;
         -0.1770 ,  0.4292  ,  0.8857 ];
    
     
     cloud2 = ( R * cloud2')';
    
     cloud1_Normal = NICP_CalculateNormal( cloud1 );
     cloud2_Normal = NICP_CalculateNormal( cloud2 );
    
     cloud1_Normal = NormalGeneral(cloud1_Normal);
     cloud2_Normal = NormalGeneral(cloud2_Normal);
     
     
     close all
     
     figure(1)
     pcshow( cloud1 , [0,1,0] , 'MarkerSize' , 40 );
     hold on
     quiver3(cloud1(:,1),cloud1(:,2),cloud1(:,3),cloud1_Normal(:,1),cloud1_Normal(:,2),cloud1_Normal(:,3));

     pcshow( cloud2 , [0,1,0] , 'MarkerSize' , 40 );
%      hold on
     quiver3(cloud2(:,1),cloud2(:,2),cloud2(:,3),cloud2_Normal(:,1),cloud2_Normal(:,2),cloud2_Normal(:,3));
     hold off
     
%      cloud2 = ( R \ cloud2')';
%      cloud2_Normal = NICP_CalculateNormal( cloud2 );
%      cloud2_Normal = NormalGeneral(cloud2_Normal);
%      
%      figure(2)
%      pcshow( cloud1 , [0,1,0] , 'MarkerSize' , 40 );
%      hold on
%      quiver3(cloud1(:,1),cloud1(:,2),cloud1(:,3),cloud1_Normal(:,1),cloud1_Normal(:,2),cloud1_Normal(:,3));
% 
%      pcshow( cloud2 , [0,1,0] , 'MarkerSize' , 40 );
% %      hold on
%      quiver3(cloud2(:,1),cloud2(:,2),cloud2(:,3),cloud2_Normal(:,1),cloud2_Normal(:,2),cloud2_Normal(:,3));
%      hold off
     
%      b = zeros(3,1);
%      
%      for i = 1 : 844
%      
%          b = b + ToAntisymmetric_Mat( cloud1_Normal(i,1:3)' ) * cloud2_Normal(i,1:3)';
%      
%      end
     
     
     NICPtestfile_normal_ortation( cloud1,cloud2,cloud1_Normal,cloud2_Normal )
     
%      RO_ICP(cloud1,cloud2)
    
     
     
     
    function [T] = RO_ICP(Refer_cloud,Match_cloud)
    
      % ���� �� ���� 1 2 �ļ��� ���أ� T �任����
    % Refer_cloud �ο����� �����任����
    % Match_cloud ƥ����� ���任����
    
    
    %% �����ʼƽ��ϵ��
    
    [ Refer_cloud_l , ~ ] = size(Refer_cloud);  %��ȡ���ƴ�С
    [ Match_cloud_l , ~ ] = size(Match_cloud);
    
    CoG_Refer_cloud = sum(Refer_cloud) / Refer_cloud_l;  %�����������
    CoG_Match_cloud = sum(Match_cloud) / Match_cloud_l;    
    
    Initial_Translation_Param = CoG_Refer_cloud - CoG_Match_cloud;  %��ʼƽ����
    Initial_Rotation_Param = diag(ones(1,3));                       %��ʼ��ת��
    
   
  %% �Ƚ��г�ʼ�任 
  
    
    
    R = Initial_Rotation_Param;   % ��ƽ���Լ���ת���г�ʼ��
    t = [0,0,0]'; 
    T(1:3,1:3) = R;
    T(1:3,4) = t;
    T(4,1:4) = [0,0,0,1];
    
    
   	%% �����ڽ���ϵ ȡŷʽ������Сֵ��Ϊ��Ӧ�� 
    
	Refer_kdtreeobj = KDTreeSearcher(Refer_cloud, 'distance' , 'euclidean' );    %���� Refercloud��kdtree
    
    
    
    
    %% ��������
    
    correspondence = zeros(Match_cloud_l,1);    
    
    
    while(true)   %�������㿪ʼ
        
    H = zeros(3);  % ���� Hessian ����
    b = zeros(3,1); % ���� J * e 
    
        for i = 1 : Match_cloud_l   %ȷ�����Ӧ��ϵ

            [ idx, ~ ] = knnsearch( Refer_kdtreeobj , Match_cloud(i,1:3) , 'dist' , 'euclidean' , 'k' , 1); 
            correspondence(i,1) = idx; 

        end
  
        
        for i = 1 : Match_cloud_l
            
            J =  - ToAntisymmetric_Mat( Match_cloud(i,1:3)' ) ; 
            
            H = H + ( J' * J ); % Hessian
            
            b = b + J' * ( Match_cloud(i,1:3)' - Refer_cloud( correspondence(i,1) , 1:3)' );
            

        end
        
        
        delt = H \ -b;         
        
        fa = delt(1:3,1);
        
        norm_delt = norm(fa);
        a = fa / norm_delt;
        
%         J_ = ( sin( norm_delt ) / norm_delt ) * eye(3) + ( 1 - sin(norm_delt) / norm_delt )*( a * a' )+(( 1 - cos(norm_delt) ) / norm_delt ) * ToAntisymmetric_Mat(a);
        deltR = cos( norm_delt ) * eye(3) + ( 1 - cos( norm_delt ) ) * ( a * a' ) + sin( norm_delt ) * ToAntisymmetric_Mat(a);
 
        Match_cloud = ( deltR * Match_cloud' )';    
        
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
    % �÷�������׼��
    
       
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
