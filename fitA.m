
x = [0;  0;  5; 10; 15; 17; 15; 20; 20; 20; 25; 30; 35; 40];
y = [0; 40;  5; 25; 10;  5; 0; 10; 25; 15;  5; 20;  8; 40];
z = [5; 10;  5; 10; 30; 30; 30;  5;  5;  5; 30;  5; 20; 10];
%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( x, y, z );

% Set up fittype and options.
ft = 'thinplateinterp';

% Fit model to data.
%[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

%{
resultMat = [0, 0, 5];
for x = 1:400
    for y = 1:400
        %if x~=1 && y~=1
            disp(x)
            x_ = x / 10;
            y_ = y / 10;
            resultMat = [resultMat; 
                         x_, y_, fitresult(x_,y_)];
        %end
    end
end
%}
resultMat=zeros(400,400);
for x = 1:400
    for y = 1:400
        disp(x)
        resultMat(x,y) = fitresult(x/10, y/10);
    end
end

%resultMat=resultMat(2:end,:);
%disp("a")
%disp(fitresult(10, 10.5))
disp(resultMat)
%plot3(resultMat(:,1), resultMat(:, 2), resultMat(:, 3))
save('fieldData', 'resultMat')
