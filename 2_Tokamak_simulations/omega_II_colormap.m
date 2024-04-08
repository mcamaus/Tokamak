%%Custom color map: IDL Omega-II colormap

T = [0,   0,   0          %// dark
     0, 0,  255         %// blue
     255, 0, 0        %// red
     255, 255, 0        %// yellow
     255, 255, 255]./255; %// white 
 %Setting color intervals length
 x = [0
     70
     130
     200
     255];
 %Interpolation between colors
 map = interp1(x/255,T,linspace(0,1,255));
 %Color bar limits form 0 to 0.7 (black to white)
 I = linspace(0,1,255);
imagesc(I(ones(1,10),:)')
colormap(map)