function [data,S]=rotation3dFourier(data,angles)

M1=[cos(angles(1)) -sin(angles(1)) 0;
    sin(angles(1))  cos(angles(1)) 0;  
    0      0  1 ;];

M2=[1 0 0;
    0 cos(angles(2)) -sin(angles(2));
    0 sin(angles(2))  cos(angles(2));];

M3=[cos(angles(3)) 0 -sin(angles(3));
    0      1  0 ;
    sin(angles(3)) 0 cos(angles(3));  
    ];

[data,S]=affineFourier(data,M3*M2*M1);
