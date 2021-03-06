colors = lines(7);
blue1 = parula(8);
blue2 = colorcube(256);
green = summer(7);
cols = [0.8*ones(1,3); 
    [blue1(1,:);[0.6 0.1 0.5];[0.5,0,0];blue2(131:131,:)] %blue2(129:131,:)];
    colors(3,:); 
    green(2:2,:); 
    [0.8 0.1 0.7];
    green(4:4,:); 
    colors(4,:);
    [0.5 0.5 0.7];
    colors(7,:)
    [0.9 0.5 0.3];
    ];

symbs = {'^','s','d'};
indication = {'melanoma','pbmc','ascites'};

save('colors.mat','cols','symbs','indication')