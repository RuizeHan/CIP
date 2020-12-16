    vec1 = data.vec_top{127+50};  % normalizaing£¿
    vec2 = data.vec_ego; 

     vec1(1,:) = vec1(1,:);
     vec2(1,:) = vec2(1,:)/(img_w/2);
     
     scale_y = data.scales_y(127+50);

     
      vec1(2,:) = vec1(2,:)./scale_y;
      vec2(2,:) = 1./vec2(2,:);
 


    figure(2)
    plot(vec1(1,:),vec1(2,:),'ro','Linewidth',1.5);
    axis([-1 1 0 max(vec1(2,:))])
    set (gcf,'Position',[950,600,600,400])



