function recessionbars(dates,ycord)
% (c) K Nimark
%introduce intervals of dates separated by commas [1982.25 1983.75 ; 1990.5 1991.0 ;]  etc
% dates=[1982.25 1983.75 ; 1990.5 1991.0 ;] ;

for j = 1:size(dates,1)
    
    x1=dates(j,1);
    x2=dates(j,2);
    y1=ycord(2);
    y2=ycord(1);
    X=[x1 x2 x2 x1 x1 ];
    Y=[y1 y1 y2 y2 y1 ];
    
    patch(X,Y,[0.9 0.9 0.9]);
    
end
    

