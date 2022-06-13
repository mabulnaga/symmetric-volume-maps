function [c_info_1, c_info_2] = create_landmarks(M1,M2)
%CREATE_LANDMARKS 
% Function to create landmarks on boundary of two meshes.
% Iteratively selects points on each surface.
% Selects first figure, then second. After selecting a point in a figure,
% press 'enter' to continue. If wanting to quit, press 'q'.

% plot the meshes
close all;
f1 = figure('units','normalized','outerposition',[0 0 0.4 1]);
trisurf(triangulation(M1.boundary_faces,M1.verts))
title('mesh 1');
axis equal;
xlabel('x');ylabel('y');zlabel('z')
f2 = figure('units','normalized','outerposition',[0.5 0 0.5 1]);
trisurf(triangulation(M2.boundary_faces,M2.verts))
title('mesh 2');
xlabel('x');ylabel('y');zlabel('z')
axis equal;
exit = false;
ct = 0;
while ~exit
   % alternate between selecting a key on each figure. after one key is pressed, then move on.
   ct = ct +1;
   for i = 1 : 2
       if(i == 1)
            figure(f1);
       else
            figure(f2);
       end
       %type 'q' to quit
%        shg
        dcm_obj = datacursormode(i);
%         set(dcm_obj,'DisplayStyle','window',...
%             'SnapToDataVertex','on','Enable','on')
       pause;
%        k = 0;
%        while(k~=1)
%            k = waitforbuttonpress;
%        end
       if strcmp(f1.CurrentCharacter ,'q') || strcmp(f2.CurrentCharacter,'q')
           exit = true;
           break;
       end
       if(i == 1)
            c_info_1(ct,:) = getCursorInfo(dcm_obj).Position;
       else
            c_info_2(ct,:) = getCursorInfo(dcm_obj).Position;
       end
   end
end


    
    
    


