function [AggregationScores,Plots,GroupingData,RootMeanSquares] = AggregationAnalysis(Uint8ImageData,AggregationSquareSize)
% Uint8ImageData: Image Data, assumed to be in Uint8 format (may work for
% other formats as well, but not tested
% AggregationSquareSize: Size of the square to select aggregation data from
% in pixels, integer

DoubleImageData = double(Uint8ImageData);
Delta = DoubleImageData(2:end,:)-DoubleImageData(1:end-1,:); %dI / dx
Delta = Delta(:,2:end)-Delta(:,1:end-1); % dI^2 / dx dy
SquareData = zeros(size(Delta,1),size(Delta,2),(AggregationSquareSize*2+1)^2);
I = 1;
for Range1 = -AggregationSquareSize:AggregationSquareSize
    for Range2 = -AggregationSquareSize:AggregationSquareSize
        SquareData(:,:,I) = circshift(Delta,[Range1,Range2]);
        I = I+1;
    end
end
DeviationData = std(SquareData,0,3);

%%%%%%%% Score Calculation

xMeanDeviation = mean(DeviationData,1);
yMeanDeviation = mean(DeviationData,2);

xMeanDeviation = xMeanDeviation-min(xMeanDeviation(:));
xMeanDeviation = xMeanDeviation./max(xMeanDeviation);
yMeanDeviation = yMeanDeviation-min(yMeanDeviation);
yMeanDeviation = yMeanDeviation./max(yMeanDeviation);

GroupingData{1} = DeviationData;
GroupingData{2} = xMeanDeviation;
GroupingData{3} = yMeanDeviation;
AggregationScores(1) = 1./sum(diff([0,find(xMeanDeviation < mean(xMeanDeviation))])>1);
AggregationScores(2) = 1./sum(diff([0;((find(yMeanDeviation < mean(yMeanDeviation))))])>1);
RootMeanSquares(1) = sqrt(sum(xMeanDeviation.^2)./numel(xMeanDeviation));
RootMeanSquares(2) = sqrt(sum(yMeanDeviation.^2)./numel(yMeanDeviation));

%%%%%%% Plot Generation
subplot(1,3,1)
hold off
Plots(1) = imagesc(Uint8ImageData);
hold on
Plots(2) = plot3(1:numel(xMeanDeviation),zeros(size(xMeanDeviation)),xMeanDeviation,'LineWidth',2);
Plots(3) = plot3(zeros(size(yMeanDeviation))+numel(xMeanDeviation),1:numel(yMeanDeviation),yMeanDeviation,'LineWidth',2);
xlabel('X')
ylabel('Y')
axis([0,numel(xMeanDeviation),0,numel(yMeanDeviation),0,max(max(yMeanDeviation),max(xMeanDeviation))])
view(3)

subplot(1,3,2)
hold off
Plots(4) = plot(xMeanDeviation);
hold on
Plots(5) = plot(xMeanDeviation.*0+mean(xMeanDeviation));
Plots(6) = plot(xMeanDeviation.*0+RootMeanSquares(1));
xlabel('X')
ylabel('Relative Intensity')
title(sprintf('X-Aggregation Score = %g',AggregationScores(1)))
axis([0,numel(xMeanDeviation),0,max(xMeanDeviation)])

subplot(1,3,3)
hold off
Plots(7) = plot(yMeanDeviation);
hold on
Plots(8) = plot(yMeanDeviation.*0+mean(yMeanDeviation));
Plots(9) = plot(yMeanDeviation.*0+RootMeanSquares(2));
xlabel('Y')
ylabel('Relative Intensity')
title(sprintf('Y-Aggregation Score = %g',AggregationScores(2)))
axis([0,numel(yMeanDeviation),0,max(yMeanDeviation)])

end