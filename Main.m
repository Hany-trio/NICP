%% 主函数 main
    %20201206 NICP 考虑法线约束的ICP算法
    clc
    clear
    
%% 读取点云数据

    cloud1 = pcread("rabbit_Segment_202011142150.pcd");
    cloud2 = pcread("rabbit_Segment_202011142150_onlyRO.pcd");   
    
    cloud1 = cloud1.Location;
    cloud2 = cloud2.Location;

    
%% 计算法线

     cloud1_Normal = NICP_CalculateNormal( cloud1 );
     cloud2_Normal = NICP_CalculateNormal( cloud2 );
    
     cloud1_Normal = NormalGeneral(cloud1_Normal);
     cloud2_Normal = NormalGeneral(cloud2_Normal);
     
	 figure(1)
     pcshow( cloud1 , [0,1,0] , 'MarkerSize' , 40 );
     hold on
     quiver3(cloud1(:,1),cloud1(:,2),cloud1(:,3),cloud1_Normal(:,1),cloud1_Normal(:,2),cloud1_Normal(:,3));
     hold off
     
     figure(2)
     pcshow( cloud2 , [0,1,0] , 'MarkerSize' , 40 );
     hold on
     quiver3(cloud2(:,1),cloud2(:,2),cloud2(:,3),cloud2_Normal(:,1),cloud2_Normal(:,2),cloud2_Normal(:,3));
     hold off
     
     close all
%% NICP
    
    HCode_NICP(cloud1,cloud2,cloud1_Normal,cloud2_Normal);
    
    
    
    
    
    
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
        norm(outNormal(i,1:3))
        
    end
    
    
    
    end
    