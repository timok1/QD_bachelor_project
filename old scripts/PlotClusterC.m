 %  29 Mar 2010
%  placing atoms in Si NC layer by layer (anion/cation/anion)
%  determining number of surface atoms
%  adding H atoms for saturation of bonds (one atom per each bond)!
%  
%  output: 3D figure
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  clc;
clear all;
opengl neverselect;


flag_carbon=1;

%  %  %  %  %  %  %  %  %  %  %  %  %  %  
  for D=[1.4];	%diameter of nanocrysral, nm
a=0.54;	%lattice constant, nm
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
path_head='data/';

D_name=sprintf('D%g_%gC%g',floor(D),floor(10.01*(D-floor(D))),flag_carbon);
path_head=[path_head,D_name,'/'];
mkdir(path_head);
points_name=[path_head,'points_',D_name,'.txt'];
cells_name=[path_head,'cells_',D_name,'.txt'];
figure_name=[path_head,'NC_3D_',D_name,'.pdf'];
disp('3D hydrogenated Si NC')
disp(['atomic coordinates will be stored in "',points_name,'"']);
disp(['elementary cell coordinates will be stored in "',cells_name,'"']);

%  %  %  %  %  %  %  %  %  %  %  %  %  %  

view_NC=[3,4,5];%view for the figure;

fontsize=17;	%size	 of font for plots
MarkerSize=30*5/(D+1);		%size of point in pixel
shift_amplitude=0.1;	%shift between two points



%  basis vectors of diamond lattice (divided by latt constant)ы
A1=1/2*[0,1,1];
A2=1/2*[1,0,1];
A3=1/2*[1,1,0];
%  nontrivial translation
shift=1/4*[1,1,1];











flag_sphere=0;	%flag specifies whether to draw boundary sphere

Amax=6;	%upper limit for cycle over the lattice sites
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  end of parameters definition
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  


points			=[];	%array to store coordinates of atoms
% x
% y
% z
% sublattice (0 or 1)
% element    (1 or 2) 
silicon  = 1;
hydrogen = 2;
carbon	=  3;
% surface    (1 or 0) - for Si only

draw_bonds		=1;
export_figure	=1;


color_C='r';
color_H='k';
color_Si='r';
r_Si	=0.15;
r_C	=0.15;	
r_H	=0.05;
n_sphere=8;

bond_color='k';bond_line_width=2;

% nearest neighbors,a/4 units
NN1=[[1,     1,    1];...
    [ 1,    -1,   -1];...
    [-1,     1,   -1];...
    [-1,    -1,    1]];
  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
%  starting adding layers after layers
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
points			=[0,0,0,0,silicon,0];
previous_layer		=points;
sublattice_previous=0;
will_add_more_layers=1;
while  will_add_more_layers
	NN=NN1*(1-2*sublattice_previous);	%nearest neigbors for this sublattice
	new_layer=[];
	for j=1:size(previous_layer,1)
	A=previous_layer(j,1:3);
	for k=1:4	% 
		  possible_neigbor = A + NN(k,:) *1/4;		
		  flag_already=false;		
		  j1=1;
		  while j1<=size(points,1)&&(~flag_already)
			if round(4*(points(j1,1:3)-possible_neigbor))==0
			  flag_already=1;end;
			j1=j1+1;
		  end
		  j1=1;
		  while j1<=size(new_layer,1)&&(~flag_already)
			if round(4*(new_layer(j1,1:3)-possible_neigbor))==0
			  flag_already=1;end;
			j1=j1+1;
		  end			
		  if ~flag_already %this atom has not been added before
			new_layer=[new_layer;possible_neigbor];
		  end
	end%possible neighbor  
  end %previous layer
  
  max_distance=0;
  siz=size(new_layer,1);
  %calculating distance between  new layer and center of NC
  for j1=1:siz
	distance=sqrt(sum(new_layer(j1,:).*new_layer(j1,:)));
	if  distance>max_distance
	  max_distance=distance;
	end
  end    
  will_add_more_layers=0;
  if max_distance*a<=D/2+sqrt(3)/4*a
	%new layer is inside NC, so it consists of Si atoms
	will_add_more_layers=1;
	points=[points;[new_layer,ones(siz,1)*(1-sublattice_previous),ones(siz,1)*silicon,zeros(siz,1)]];
%  	for j1=1:siz 	 colors=[colors,c_Si];end
	
  else %new layer is outside
	%new layer is outside NC, so it consists of hydrogen atoms
	will_add_more_layers=0;
	for j=1:size(previous_layer,1)
	  A=previous_layer(j,1:3);
	  for k=1:4	% 
			possible_H = A + NN(k,:) *1/4;	
			possible_H_to_put = A + NN(k,:) *1/4*1/3;		
			flag_H=false;		
			j1=1;
			while j1<=size(new_layer,1)&&(~flag_H)
			  if round(4*(new_layer(j1,1:3)-possible_H))==0
				flag_H=1;end;
				j1=j1+1;
			end
			if flag_H
			 points=[points;[possible_H_to_put,1-sublattice_previous,hydrogen,0]];
%  			  colors=[colors,c_H];
			end
	  end
	end
%  %  	  points=[points;[new_layer,ones(siz,1)*(1-sublattice_previous),ones(siz,1)*hydrogen,zeros(siz,1)]];
%  %  	  for j1=1:siz 	 colors=[colors,c_H];end

%  NB: here we need one H atom per each bond! so the number of H atoms is in fact larger than number of Si atoms???



%  	  %previous layer of Si atoms is marked as boundary	  
	  
	  points(previous_indices,6)=1;
%  	  colors(previous_indices)=c_Si_boundary;	  
	
  end  %if new layer is outside
  previous_indices=[size(points,1)-siz+1:size(points,1)];
  sublattice_previous=1-sublattice_previous;
  previous_layer=new_layer;
  
end %iteration to add layer of atoms

NNCA=size(points,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  adding carbon
if flag_carbon>0
  atoms_Si=find(points(:,5)==silicon);
  [m,ind]=max(points(atoms_Si,1).^2+points(atoms_Si,2).^2+points(atoms_Si,3).^2);
  points(ind,5)=carbon;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%starting plotting results

%ind=find(points(:,5)==hydrogen);
%points=points(ind,:);
%colors=colors(ind);
%NNCA=numel(ind);

figure(1);clf(gcf);
N_Si=0;
N_H=0;
N_Si_bound=0;
N_C=0;
for i=1:NNCA
    %for i=1:100
%      scatter3(a*points(i,1),a*points(i,2),a*points(i,3),'filled','MarkerFaceColor',colors(i),'SizeData',MarkerSize);		hold on;


%  if points(j,4)==Pb
	  
%  	else
%  	  r=r_Se;color=color_Se;
%  	end
 	

    if points(i,5) == silicon
        N_Si=N_Si+1;
        r=r_Si;color=color_Si;
    end
    if points(i,5) == carbon
        N_C=N_C+1;
        r=r_C;color=color_C;
    end
    if points(i,5) == hydrogen
        N_H=N_H+1;
        r=r_H;color=color_H;
    end
    
    if points(i,6) == 1
        N_Si_bound=N_Si_bound+1;
    end
    
    [x,y,z]=ellipsoid(points(i,1),points(i,2),points(i,3),r,r,r,n_sphere);
    p=surf(x,y,z);
    set(p,'FaceColor',color,'EdgeColor','none');
    hold on;
end % for

if draw_bonds
%  	2+2
%  	  are_neigbors=sparse(zeros(N_atoms));
	  for k=1:NNCA
%  	  keyboard
		for l=k+1:NNCA
		  r_k=points(k,1:3);
		  r_l=points(l,1:3);		  
		  
%  		  if sum(abs(r_k-r_l))==1
  		  if abs(sum(abs(r_k-r_l).^2)-3/16)<1e-4...
  		  &&points(k,5)~=hydrogen&&points(l,5)~=hydrogen
%    		  keyboard;
			p=line([r_k(1),r_l(1)],[r_k(2),r_l(2)],[r_k(3),r_l(3)]);
			set(p,'LineWidth',bond_line_width,'Color',bond_color);
		  end
		end
	  end	  
end %draw_bonds

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
 if flag_sphere % 1st sphere 	
    [X,Y,Z] =sphere;
    X=X*D/(2);Y=Y*D/(2);Z=Z*D/(2);
      mesh(X,Y,Z,'FaceColor','none','EdgeColor',c_Si_boundary);
      hold on;
  end
set(gca,'fontsize',fontsize);
%  axis off;hold off;
axis tight;daspect([1 1 1]);
   axis off;
     camlight(320,50); lighting phong;
     




xlabel('x, nm','fontsize',fontsize);
ylabel('y, nm','fontsize',fontsize);
zlabel('z, nm','fontsize',fontsize);
title([sprintf('D=%g nm',D), ' , N_{atom}^{(tot)}=',sprintf('%g',NNCA),', \newline  N_{Si}^{(tot)}=',sprintf('%g',N_Si),' , N_{H}^{(tot)}=',sprintf('%g',N_H),' , N_{Si}^{(surf)}=',sprintf('%g',N_Si_bound)],'fontsize',fontsize);
disp([sprintf('D=%g nm',D), ' , N_{atom}^{(tot)}=',sprintf('%g',NNCA),', N_{Si}^{(tot)}=',sprintf('%g',N_Si),' , N_{H}^{(tot)}=',sprintf('%g',N_H),' , N_{Si}^{(surf)}=',sprintf('%g',N_Si_bound)])
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %   
view(view_NC);
dlmwrite(points_name,points);
%  dlmwrite(cells_name,cells);
print('-dpdf',figure_name);

end %Diameter